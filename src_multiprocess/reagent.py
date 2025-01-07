import numpy as np
from rdkit import Chem
import math


class Reagent:
    def __init__(self, reagent_name: str, smiles: str):
        """
        Basic init
        :param reagent_name: Reagent name
        :param smiles: smiles string
        """
        self.reagent_name = reagent_name
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(self.smiles)
        self.scores_batch = []
        # Will be initialized during init_given_prior
        self.var_score_known = None  
        self.std_score_known = None
        self.posterior_mean = None
        self.posterior_std = None
        self.sum_w = None

    def init_prior(self, prior_mean_score: float, prior_std_score: float):
        """
        From the score distribution seen across all reagents during the warm up phase.
        """
        self.var_score_known = prior_std_score ** 2
        self.posterior_std = prior_std_score
        self.posterior_mean = prior_mean_score
        self.std_score_known = prior_std_score 
        self.sum_w = math.exp(self.posterior_mean/self.std_score_known)

    def single_update(self, observed_value: float):
        """
        Does the bayesian update of the posterior mean and standard deviation.
        :param observed_value: New score collected for the reagent
        """
        # The posterior variance now serve as the prior
        prior_var = self.posterior_std ** 2
        denominator = prior_var + self.var_score_known
        # mean --> Boltzmann weighted average
        w = math.exp(observed_value/self.std_score_known)
        self.sum_w += w
        self.posterior_mean = self.posterior_mean + (w/self.sum_w)*(observed_value-self.posterior_mean)
        # std
        numerator = prior_var * self.var_score_known
        self.posterior_std = np.sqrt(numerator / denominator)
        
    def add_score(self, score: float):
        self.scores_batch.append(score)

    def multiple_update(self):
        if self.scores_batch:
           n = len(self.scores_batch)
           prior_var = self.posterior_std ** 2
           denominator = n * prior_var + self.var_score_known
           # mean --> Boltzmann weighted average
           scores_batch = np.array(self.scores_batch)
           w_batch = np.exp(scores_batch/self.std_score_known)
           mean_batch = np.average(scores_batch,weights=w_batch)
           w_sum_batch = np.sum(w_batch)
           self.sum_w += w_sum_batch
           self.posterior_mean = self.posterior_mean + (w_sum_batch/self.sum_w)*(mean_batch-self.posterior_mean)
           # std
           numerator = prior_var * self.var_score_known
           self.posterior_std = np.sqrt(numerator / denominator)   
           # reset scores_batch for next batch update    
           self.scores_batch = []
            
