# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 16:13:35 2023

@author: 20213628
"""

from rdkit import Chem
from rdkit.Contrib.SA_Score import sascorer
import seaborn as sns
import matplotlib.pyplot as plt

def list_sa_scores_distribution_with_canonicalization(smiles_list):
    """
    Plot the distribution of Synthetic Accessibility (SA) scores for a list of molecules after canonicalization.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Output:
    - sa_scores (list) : List of sa_scores
    """

    sa_scores = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Canonicalize SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            # Calculate SA score using the sascorer function
            sa_score = sascorer.calculateScore(Chem.MolFromSmiles(canonical_smiles))
            sa_scores.append(sa_score)

    return sa_scores



