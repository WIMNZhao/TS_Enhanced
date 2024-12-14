import random
from multiprocessing import Pool
from typing import List, Optional, Tuple
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm.auto import tqdm
from reagent import Reagent
from ts_logger import get_logger
from ts_utils import read_reagents
from evaluators import DBEvaluator


class ThompsonSampler:
    def __init__(self, processes: int, scaling: float, log_filename: Optional[str] = None):
        """
        Basic init
        :param log_filename: Optional filename to write logging to.
        """
        # A list of lists of Reagents. Each component in the reaction will have one list of Reagents in this list
        self.reagent_lists: List[List[Reagent]] = []
        self.reaction = None
        self.evaluator = None
        self.num_prods = 0
        self.processes = processes
        self.scaling = scaling
        self.hide_progress = False
        self.num_warmup = None
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
        return [score, product_smiles, product_name]

    def warm_up(self, num_warmup_trials=3, results_filename="results.csv"):
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
            if np.isfinite(r[0]):             
               warmup_scores.append(r[0])
               for idx, rdx in enumerate(p):
                   self.reagent_lists[idx][rdx].add_score(r[0]*self.scaling)
        # initialize each reagent
        prior_mean = np.mean(warmup_scores)
        prior_std = np.std(warmup_scores)
        for i in range(0, len(self.reagent_lists)):
            for j in range(0, len(self.reagent_lists[i])):
                reagent = self.reagent_lists[i][j]
                reagent.init_prior(prior_mean_score=prior_mean, prior_std_score=prior_std)
                reagent.multiple_update()

        # write the results to disk
        out_df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])
        out_df.to_csv(results_filename, index=False, na_rep='nan')
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
        self.num_warmup = len(results)
        return self.num_warmup

    def search(self, num_per_cycle:int=1000, percent_of_library:float=0.001, stop:int=1000,results_filename="results.csv"):
        """
        Run the search with roulette wheel selection 
        :param: scaling: scale the temperature; positive if maximum score is preferred negative otherwise
        :param: num_per_cycle: the number of products to make per search cycle
        :param: percent_of_library: percent of the library to be screened
        :param: stop: number of resamples on a roll
        :return: a list of SMILES and scores
        """
        uniq = {}
        out_list = []
        n_resample = 0
        rng = np.random.default_rng()
        nsearch = int(percent_of_library*self.num_prods - self.num_warmup)
        count = 0
        
        idxs_component = range(0,len(self.reagent_lists))

        with tqdm(total=nsearch, bar_format='{l_bar}{bar}| {elapsed}<{remaining}, {rate_fmt}{postfix}', disable=self.hide_progress) as pbar:
          while (len(uniq) < nsearch):
            matrix = []
            pairs = []

            # thermal cycling between normal and greedy-selection adapted roulette wheel selection
            if random.uniform(0, 1) < 0.8:   
               ttt = random.uniform(7,9)
               idx_c = random.choice(idxs_component)
               app_tc = True
            else:
               app_tc = False

            for ii, rg in enumerate(self.reagent_lists):
                # sample scores
                stds = np.array([r.posterior_std  for r in rg])
                mu   = np.array([r.posterior_mean for r in rg])
                rg_score = rng.normal(size=len(rg)) * stds + mu
                # apply thermal cycling; room temp is std of the sampled scores per reaction component
                if app_tc:
                   # increase temp for one component 
                   if ii == idx_c:
                      rg_score = np.exp(rg_score/(np.std(rg_score)*ttt))
                   # decrease temp for others   
                   else:
                      rg_score = np.exp(rg_score/np.std(rg_score)*ttt)
                else:   
                   rg_score = np.exp(rg_score/np.std(rg_score))
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
                if np.isfinite(r[0]):
                   out_list.append(r)
                   for idx, rdx in enumerate(p):
                       self.reagent_lists[idx][rdx].single_update(r[0]*self.scaling)                          
            # write to disk
            out_df = pd.DataFrame(results)
            out_df.to_csv(results_filename, mode='a', header=False, index=False, na_rep='nan')
            # stop criteria check
            if n_resample >= stop: 
                self.logger.info(f"Stop criteria reached with the number of resamples: {n_resample}")
                break
            # logging
            if count % 100 == 0:
                if self.scaling > 0:
                   top_score, top_smiles, top_name = max(out_list)
                else:
                   top_score, top_smiles, top_name = min(out_list)
                self.logger.info(f"Iteration: {count} best score: {top_score:2f} smiles: {top_smiles} combi: {top_name}")
            # progressbar update
            incr = len(pairs_u)
            if pbar.n + incr > nsearch:
               incr = nsearch - pbar.n
            pbar.update(incr)
            # iterations
            count += 1
        return out_list

