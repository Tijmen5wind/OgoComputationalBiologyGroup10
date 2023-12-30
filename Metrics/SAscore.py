# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 16:13:35 2023

@author: 20213628
"""

from rdkit import Chem
from rdkit.Contrib.SA_Score import sascorer

def calculate_average_sa_score_with_canonicalization(smiles_list):
    """
    Calculate the average Synthetic Accessibility (SA) score for a list of molecules after canonicalization.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - average_sa_score (float): Average Synthetic Accessibility score.
    """

    total_sa_score = 0.0
    valid_count = 0

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Canonicalize SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            # Calculate SA score using the sascorer function
            sa_score = sascorer.calculateScore(Chem.MolFromSmiles(canonical_smiles))
            total_sa_score += sa_score
            valid_count += 1

    if valid_count == 0:
        return 0.0

    average_sa_score = total_sa_score / valid_count
    return average_sa_score

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "CCO", "C1=CC=CC=C1", "HOI"]
average_sa_score = calculate_average_sa_score_with_canonicalization(smiles_list_to_check)

print(f"Average SA score: {average_sa_score:.2f}")

