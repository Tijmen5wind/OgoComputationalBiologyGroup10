# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 23:35:15 2023

@author: 20213628
"""

import pandas as pd
from rdkit import Chem

def validity_percentage(smiles_list):
    """
    Calculate the percentage of valid SMILES strings.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - percentage_valid (float): Percentage of valid SMILES strings.
    """

    df = pd.DataFrame({"SMILES": smiles_list})
    df["Mol"] = df["SMILES"].apply(Chem.MolFromSmiles)

    total_smiles = len(df)
    valid_smiles = df["Mol"].count()

    if total_smiles == 0:
        return 0.0

    percentage_valid =valid_smiles / 1000
    return percentage_valid

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "InvalidSMILES", "C1CC1"]

percentage_valid = validity_percentage(smiles_list_to_check)

print(f"Validity Percentage: {percentage_valid:.2f}%")


