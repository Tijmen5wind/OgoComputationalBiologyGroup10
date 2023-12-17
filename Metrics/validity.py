# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 23:35:15 2023

@author: 20213628
"""

from rdkit import Chem

def validity_percentage(smiles_list):
    """
    Calculate the percentage of valid SMILES strings.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - percentage_valid (float): Percentage of valid SMILES strings.
    """

    total_smiles = len(smiles_list)
    valid_smiles = 0

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            valid_smiles += 1

    if total_smiles == 0:
        return 0.0

    percentage_valid = (valid_smiles / total_smiles) * 100
    return percentage_valid

# Example usage:
smiles_to_check = ["CCO", "C1=CC=CC=C1", "CC", "CC(C)(C)C(=O)O"]

percentage_valid = validity_percentage(smiles_to_check)

print(f"Percentage of valid SMILES: {percentage_valid:.2f}%")

