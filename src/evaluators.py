from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

class FPEvaluator:
    """An evaluator class that calculates a fingerprint Tanimoto to a reference molecule
    """
    def __init__(self, input_dict):
        self.ref_smiles = input_dict["query_smiles"]
        self.ref_fp = self.mol2morgan_fp(Chem.MolFromSmiles(self.ref_smiles))
        self.num_evaluations = 0

    @property
    def counter(self):
        return self.num_evaluations

    def mol2morgan_fp(self, mol):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        return fp

    def evaluate(self, rd_mol_in):
        self.num_evaluations += 1
        rd_mol_fp = self.mol2morgan_fp(rd_mol_in)
        return DataStructs.TanimotoSimilarity(self.ref_fp, rd_mol_fp)


