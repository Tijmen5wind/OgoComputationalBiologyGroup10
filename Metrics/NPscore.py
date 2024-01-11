# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 16:27:55 2023

@author: 20213628
"""
from rdkit import Chem
from rdkit.Contrib.NP_Score import npscorer
import seaborn as sns
import matplotlib.pyplot as plt

def list_np_scores_distribution_with_canonicalization(smiles_list):
    """
    List the distribution of Natural Product (NP) scores for a list of molecules after canonicalization.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Output:
    - np_scores (list): List of NP scores
    """

    np_scores = []
    fscore = npscorer.readNPModel()

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Canonicalize SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            # Calculate NP score using the npscorer function
            np_score, confidence = npscorer.scoreMolWConfidence(Chem.MolFromSmiles(canonical_smiles), fscore)
            np_scores.append(np_score)
    return np_scores










