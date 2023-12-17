# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 00:25:29 2023

@author: 20213628
"""

from rdkit import Chem
from rdkit.Chem import inchi

def calculate_novelty(smiles_list, dataset):
    """
    Calculate the novelty of SMILES strings based on a dataset (chembl)

    Parameters:
    - smiles_list (list): List of SMILES strings to check for novelty.
    - dataset (list): Reference dataset of SMILES strings.

    Returns:
    - novelty_percentage (float): Percentage of novel SMILES strings.
    """

    unique_smiles = set(smiles_list)
    unique_reference = set(dataset)

    novel_smiles = unique_smiles - unique_reference

    if not unique_smiles:
        return 0.0

    novelty_percentage = (len(novel_smiles) / len(unique_smiles)) * 100
    return novelty_percentage

# Example usage:
reference_dataset = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O"]

smiles_to_check = ["CCO", "C1=CC=CC=C1", "InvalidSMILES", "CC(C)(C)C(=O)O", "CCN"]

novelty_percentage = calculate_novelty(smiles_to_check, reference_dataset)

print(f"Novelty Percentage: {novelty_percentage:.2f}%")
