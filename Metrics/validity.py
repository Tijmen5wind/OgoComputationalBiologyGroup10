# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 23:35:15 2023

@author: 20213628
"""

from rdkit import Chem

def calculate_validity(smiles_list):
    """
    Calculate the validity of SMILES strings.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - validity (float): validity of the SMILES strings.
    """

    total_smiles = len(smiles_list)
    valid_smiles = 0

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            valid_smiles += 1

    if total_smiles == 0:
        return 0.0

    validity = (valid_smiles / total_smiles) 
    return validity

# Example usage:
smiles_to_check = ["CCO", "C1=CC=CC=C1", "CC", "CC(C)(C)C(=O)O"]

validity = calculate_validity(smiles_to_check)

print(validity)

