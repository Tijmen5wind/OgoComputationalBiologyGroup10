# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 15:53:31 2023

@author: 20213628
"""

def uniqueness(smiles_list):
    """
    Count the percentage of unique SMILES in a list.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - uniqueness (float): percentage of unique SMILES.
    """

    unique_smiles_set = set(smiles_list)
    unique_count = len(unique_smiles_set)
    uniqueness=unique_count/len(smiles_list)

    return uniqueness

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "CCO", "C1=CC=CC=C1"]

uniqueness = uniqueness(smiles_list_to_check)

print(uniqueness)
