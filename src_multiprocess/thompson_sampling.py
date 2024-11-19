import random, math
from multiprocessing import Pool
from typing import List, Optional, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm.auto import tqdm
from reagent import Reagent
from ts_logger import get_logger
from ts_utils import read_reagents
from evaluators import DBEvaluator


class ThompsonSampler:
    def __init__(self, processes: int, log_filename: Optional[str] = None):
        """
        Basic init
        :param log_filename: Optional filename to write logging to.
        """
        # A list of lists of Reagents. Each component in the reaction will have one list of Reagents in this list
        self.reagent_lists: List[List[Reagent]] = []
        self.reaction = None
        self.evaluator = None
        self.num_prods = 0
        self.T = None
        self.processes = processes
        self.hide_progress = False
        #self.fail_score = 0
        self.logger = get_logger(__name__, filename=log_filename)

    def set_hide_progress(self, hide_progress: bool) -> None:
        """
        Hide the progress bars
        :param hide_progress: set to True to hide the progress baars
        """
        self.hide_progress = hide_progress

    def read_reagents(self, reagent_file_list, num_to_select: Optional[int] = None):
        """
        Reads the reagents from reagent_file_list
        :param reagent_file_list: List of reagent filepaths
        :param num_to_select: Max number of reagents to select from the reagents file (for dev purposes only)
        :return: None
        """
        self.reagent_lists = read_reagents(reagent_file_list, num_to_select)
        self.num_prods = np.prod([len(x) for x in self.reagent_lists])
        self.logger.info(f"{self.num_prods:.2e} possible products")
        
    def get_num_prods(self) -> int:
        """
        Get the total number of possible products
        :return: num_prods
        """
        return self.num_prods

    def set_evaluator(self, evaluator):
        """
        Define the evaluator
        :param evaluator: evaluator class, must define an evaluate method that takes an RDKit molecule
        """
        self.evaluator = evaluator

    def set_reaction(self, rxn_smarts):
        """
        Define the reaction
        :param rxn_smarts: reaction SMARTS
        """
        self.reaction = AllChem.ReactionFromSmarts(rxn_smarts)

    def evaluate(self, choice_list: List[int]) -> Tuple[str, float]:
        """Evaluate a set of reagents
        :param choice_list: list of reagent ids
        :return: smiles for the reaction product, score for the reaction product
        """
        selected_reagents = []
        for idx, choice in enumerate(choice_list):
            selected_reagents.append(self.reagent_lists[idx][choice])
        prod = self.reaction.RunReactants([reagent.mol for reagent in selected_reagents])
        product_name = "_".join([reagent.reagent_name for reagent in selected_reagents])
        score = np.nan
        product_smiles = "FAIL"
        if prod:
           prod_mol = prod[0][0]  # RunReactants returns Tuple[Tuple[Mol]]
           Chem.SanitizeMol(prod_mol)
           product_smiles = Chem.MolToSmiles(prod_mol)
           if isinstance(self.evaluator, DBEvaluator):
              score = self.evaluator.evaluate(product_name)
              score = float(score)
           else:
              score = self.evaluator.evaluate(prod_mol)
        return [product_smiles, product_name, score]

    def warm_up(self, num_warmup_trials=3):
        """
        Warm-up phase by stochastic paralell pairing repeated for num_warmup_trials
        :param num_warmup_trials: number of warm-up cycles
        """
        pairs = []
        reagent_count_list = [len(x) for x in self.reagent_lists]
        idx = np.argmax(reagent_count_list)
        rmax = reagent_count_list[idx]

        pack = []
        for i in reagent_count_list:
            if i < rmax:
                pack.append((1,rmax // i,rmax % i))
            else:
                pack.append((0,))    

        for i in range(num_warmup_trials):
            matrix = []
            for j, nr in enumerate(reagent_count_list):
                idx_r = list(range(0,nr))
                random.shuffle(idx_r)
                if pack[j][0]:
                   matrix.append(idx_r * pack[j][1] + idx_r[:pack[j][2]] )
                else:
                   matrix.append(idx_r)
            pairs.extend(np.array(matrix).transpose())  
        #
        warmup_scores = []
        if self.processes >1:
           with Pool(self.processes) as pool:
                results = pool.map(self.evaluate, pairs)
        else:
           results = []
           for p in pairs:
               results.append(self.evaluate(p))
        # avoid multiprocessing overhead if nprocess == 1 for fast scoring method
        for r, p in zip(results,pairs):
            if np.isfinite(r[2]):             
               warmup_scores.append(r[2])
               for idx, rdx in enumerate(p):
                   self.reagent_lists[idx][rdx].add_score(r[2])
        # initialize each reagent
        prior_mean = np.mean(warmup_scores)
        prior_std = np.std(warmup_scores)
        for i in range(0, len(self.reagent_lists)):
            for j in range(0, len(self.reagent_lists[i])):
                reagent = self.reagent_lists[i][j]
                reagent.init_prior(prior_mean_score=prior_mean, prior_std_score=prior_std)
                reagent.multiple_update()
        # logging
        self.logger.info(
            f"warmup score stats: "
            f"cnt={len(warmup_scores)}, "
            f"mean={np.mean(warmup_scores):0.4f}, "
            f"std={np.std(warmup_scores):0.4f}, "
            f"min={np.min(warmup_scores):0.4f}, "
            f"max={np.max(warmup_scores):0.4f}"
            )
        # set Temperature
        self.T = prior_std

    def search(self, scaling:float=1, ld:float=0.1, decay:float=1, dcycles:int=100, num_per_cycle:int=100, num_ts_iterations:int=1000, stop:int=100):
        """
        Run the search with roulette wheel selection 
        :param: scaling: scale the temperature; positive if maximum score is preferred negative otherwise
        :param: decay: decrease temperature every dcycles by mulitiplying the current T
        :param: dcycles: the number of cycles to decrease temperature
        :param: num_per_cycle: the number of products to make per search cycle
        :param: num_ts_iterations: number of search iterations
        :param: stop: number of resamples on a roll
        :param: ld: lower bound of the temperature
        :return: a list of SMILES and scores
        """
        Temp = self.T * scaling
        uniq = {}
        out_list = []
        n_resample = 0
        ld_reached = False

        rng = np.random.default_rng()
        for i in tqdm(range(0, num_ts_iterations), desc="Cycle", disable=self.hide_progress):
            matrix = []
            pairs = []
            for rg in self.reagent_lists:
                # updated reagent score sampling from Pat Walters
                stds = np.array([r.posterior_std  for r in rg])
                mu   = np.array([r.posterior_mean for r in rg])
                rg_score = np.exp((rng.normal(size=len(rg)) * stds + mu)/Temp)
                # roulette wheel selection
                sele = np.random.choice(len(rg),num_per_cycle,p=rg_score / np.sum(rg_score))
                matrix.append(sele)
            pairs = np.array(matrix).transpose()

            pairs_u = []
            for comb in pairs:
                ss = '_'.join(str(rr) for rr in comb)
                if ss not in uniq:
                   pairs_u.append(comb)
                   uniq[ss] = None
                   n_resample = 0
                else:
                   n_resample += 1        
            #
            if self.processes >1:
               with Pool(self.processes) as pool:
                    results = pool.map(self.evaluate, pairs_u)
            else:
               results = []
               for p in pairs_u:
                   results.append(self.evaluate(p))
            # avoid multiprocessing overhead if nprocess == 1 for fast scoring method        
            for r, p in zip(results,pairs_u):
                if np.isfinite(r[2]):
                   out_list.append([r[2], r[0], r[1]])
                   for idx, rdx in enumerate(p):
                       self.reagent_lists[idx][rdx].single_update(r[2])                          

            # stop criteria check
            if n_resample == stop: 
                self.logger.info(f"Stop criteria reached with the number of resamples: {n_resample}")
                break

            if not ld_reached and i % dcycles == 0:
                Temp *= decay
                if Temp < self.T * ld:
                    Temp = self.T * ld
                    ld_reached = True

            if i % 100 == 0:
                if scaling > 0:
                   top_score, top_smiles, top_name = max(out_list)
                else:
                   top_score, top_smiles, top_name = min(out_list)
                self.logger.info(f"Iteration: {i} best score: {top_score:2f} smiles: {top_smiles} combi: {top_name}")
        return out_list
