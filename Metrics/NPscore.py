# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 16:27:55 2023

@author: 20213628
"""
from rdkit import Chem
from rdkit.Contrib.NP_Score import npscorer

def calculate_average_np_score_with_canonicalization(smiles_list):
    """
    Calculate the average Natural Product (NP) score for a list of molecules after canonicalization.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - average_np_score (float): Average Natural Product score.
    """

    total_np_score = 0.0
    valid_count = 0
    fscore = npscorer.readNPModel()

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Canonicalize SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            # Calculate NP score using the npscorer function
            np_score, confidence = npscorer.scoreMolWConfidence(Chem.MolFromSmiles(canonical_smiles), fscore)
            total_np_score += np_score
            valid_count += 1

    if valid_count == 0:
        return 0.0

    average_np_score = total_np_score / valid_count
    return average_np_score

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "CCO", "C1=CC=CC=C1"]
average_np_score = calculate_average_np_score_with_canonicalization(smiles_list_to_check)

print(f"Average NP score: {average_np_score:.2f}")




