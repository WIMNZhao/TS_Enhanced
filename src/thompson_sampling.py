import random
from typing import List, Optional, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm.auto import tqdm

from reagent import Reagent
from ts_logger import get_logger
from ts_utils import read_reagents


class ThompsonSampler:
    def __init__(self, log_filename: Optional[str] = None):
        """
        Basic init
        :param log_filename: Optional filename to write logging to. If None, logging will be output to stdout
        """
        # A list of lists of Reagents. Each component in the reaction will have one list of Reagents in this list
        self.reagent_lists: List[List[Reagent]] = []
        self.reaction = None
        self.evaluator = None
        self.num_prods = 0
        self.T = None
        self.logger = get_logger(__name__, filename=log_filename)
        self.hide_progress = False

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
        res = np.nan
        product_smiles = "FAIL"
        if prod:
            prod_mol = prod[0][0]  # RunReactants returns Tuple[Tuple[Mol]]
            Chem.SanitizeMol(prod_mol)
            product_smiles = Chem.MolToSmiles(prod_mol)
            res = self.evaluator.evaluate(prod_mol)
            if np.isfinite(res):
               [reagent.add_score(res) for reagent in selected_reagents]
        return product_smiles, product_name, res

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

        warmup_scores = []
        for p in pairs:
            _, _, score = self.evaluate(p)
            if np.isfinite(score):
               warmup_scores.append(score)
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

    def search(self, scaling:float=1, ld:float=0.1, decay:float=1, dcycles:int=100, num_per_cycle:int=3, num_cycles:int=1000, stop:int=100):
        """
        Run the search with roulette wheel selection 
        :param: scaling: scale the temperature; positive if maximum score is preferred negative otherwise
        :param: decay: decrease temperature every dcycles by mulitiplying the current T
        :param: dcycles: the number of cycles to decrease temperature
        :param: num_per_cycle: the number of products to make per search cycle
        :param: num_cycles: number of search iterations
        :param: stop: number of resamples on a roll
        :param: ld: lower bound of the temperature
        :return: a list of SMILES and scores
        """
        Temp = self.T * scaling
        uniq = {}
        out_list = []
        n_resample = 0
        ld_reached = False
        for i in tqdm(range(0, num_cycles), desc="Cycle", disable=self.hide_progress):
            matrix = []
            pairs = []
            for rg in self.reagent_lists:
                rg_score = np.zeros(len(rg))  # Create a list of scores for each reagent
                for reagent_idx, reagent in enumerate(rg):
                    rg_score[reagent_idx] = np.exp(reagent.sample() / Temp)
                sele = np.random.choice(len(rg),num_per_cycle,p=rg_score / np.sum(rg_score))
                matrix.append(sele)
            pairs = np.array(matrix).transpose()

            for ii, comb in enumerate(pairs):
                ss = '_'.join(str(rr) for rr in comb)
                if ss not in uniq:
                    smiles, name, score = self.evaluate(comb)
                    if np.isfinite(score):
                       for idx, rdx in enumerate(comb):
                           self.reagent_lists[idx][rdx].multiple_update()
                       sample_efficiency = (len(out_list) + 1) / (i*num_per_cycle+ii+1) * 100                            
                       out_list.append([score, smiles, name + '_%.1f'%sample_efficiency])
                       n_resample = 0
                       uniq[ss] = None
                else:
                    n_resample += 1

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
