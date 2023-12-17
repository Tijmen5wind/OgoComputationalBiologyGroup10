# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 19:00:27 2023

@author: 20213628
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

def lipinski_descriptors(smiles):
    """Calculates the descriptors that are within the Lipinski rule of 5 and returns them as a dictionary"""
    mol = Chem.MolFromSmiles(smiles)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    return {'Hydrogen Bond Donors': hbd, 'Hydrogen Bond Acceptors': hba, 'Molecular Weight': mw, 'LogP': logp}

def calculate_average_lipinski_properties(smiles_list):
    """
    Calculate the average Lipinski Rule of 5 properties for a list of SMILES strings.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - average_properties (dict): Dictionary with average Lipinski Rule of 5 properties.
    """
    properties_list = [lipinski_descriptors(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]

    if not properties_list:
        return None

    num_molecules = len(properties_list)

    # Initialize sums
    sum_dict = {prop: 0.0 for prop in properties_list[0]}

    # Sum up values
    for properties in properties_list:
        for prop, value in properties.items():
            sum_dict[prop] += value

    # Calculate averages
    average_properties = {prop: sum_value / num_molecules for prop, sum_value in sum_dict.items()}

    return average_properties

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "C(CCCCCCCC)CCCCCCCC(=O)O"]

average_lipinski_properties = calculate_average_lipinski_properties(smiles_list_to_check)

print(average_lipinski_properties)
