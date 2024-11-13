import numpy as np
from rdkit import Chem


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
        self.posterior_mean = None
        self.posterior_std = None

    def init_prior(self, prior_mean_score: float, prior_std_score: float):
        """
        From the score distribution seen across all reagents during the warm up phase.
        """
        self.var_score_known = prior_std_score ** 2
        self.posterior_mean = prior_mean_score
        self.posterior_std = prior_std_score

    def sample(self) -> float:
        """
        Takes a random sample from the posterior distribution of mean
        :return: sample from the posterior distribution
        """
        return np.random.normal(loc=self.posterior_mean, scale=self.posterior_std)

    def single_update(self, observed_value: float):
        """
        Does the bayesian update of the posterior mean and standard deviation.
        :param observed_value: New score collected for the reagent
        """
        # The posterior variance now serve as the prior
        prior_var = self.posterior_std ** 2
        denominator = prior_var + self.var_score_known
        # mean
        numerator = prior_var * observed_value + self.var_score_known * self.posterior_mean
        self.posterior_mean = numerator / denominator
        # std
        numerator = prior_var * self.var_score_known
        self.posterior_std = np.sqrt(numerator / denominator)
        
    def add_score(self, score: float):
        self.scores_batch.append(score)

    def multiple_update(self):
        n = len(self.scores_batch)
        s = np.sum(self.scores_batch)
        prior_var = self.posterior_std ** 2
        denominator = n * prior_var + self.var_score_known
        # mean
        numerator = s * prior_var + self.var_score_known * self.posterior_mean
        self.posterior_mean = numerator / denominator
        # std
        numerator = prior_var * self.var_score_known
        self.posterior_std = np.sqrt(numerator / denominator)   
        # reset scores_batch for next batch update    
        self.scores_batch = []
            