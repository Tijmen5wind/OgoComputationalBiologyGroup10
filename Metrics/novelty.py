# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 00:25:29 2023

@author: 20213628
"""


def calculate_novelty(smiles_list, dataset):
    """
    Calculate the novelty of SMILES strings based on a dataset (chembl)

    Parameters:
    - smiles_list (list): List of SMILES strings to check for novelty.
    - dataset (list): Reference dataset of SMILES strings.

    Returns:
    - novelty_percentage (float): Percentage of novel SMILES strings.
    """

    novel_smiles=0
    for smile in smiles_list:
        if smile not in dataset:
            novel_smiles+=1
            
    novelty_percentage = (novel_smiles) / len(smiles_list)
    return novelty_percentage

# Example usage:
reference_dataset = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O"]

smiles_to_check = ["CCO", "C1=CC=CC=C1", "InvalidSMILES", "CC(C)(C)C(=O)O", "CCN", "CCN","CCO"]

novelty_percentage = calculate_novelty(smiles_to_check, reference_dataset)

print(f"Novelty Percentage: {novelty_percentage:.2f}%")
