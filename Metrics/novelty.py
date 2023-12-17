# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 00:25:29 2023

@author: 20213628
"""

from rdkit import Chem

def novelty(smiles_list, dataset):
    """
    Calculate the novelty of SMILES strings based on a dataset (chembl)

    Parameters:
    - smiles_list (list): List of SMILES strings to check for novelty.
    - dataset (list): Reference dataset of SMILES strings.

    Returns:
    - novelty_fraction (float): Percentage of novel SMILES strings.
    """

    unique_smiles = set(smiles_list)
    unique_reference = set(dataset)

    novel_smiles = unique_smiles - unique_reference

    if not unique_smiles:
        return 0.0

    novelty_fraction = (len(novel_smiles) / len(unique_smiles))
    return novelty_fraction

# Example usage:
reference_dataset = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O"]

smiles_to_check = ["CCO", "C1=CC=CC=C1", "InvalidSMILES", "CC(C)(C)C(=O)O", "CCN"]

novelty_fraction = novelty(smiles_to_check, reference_dataset)

print(novelty_fraction)
