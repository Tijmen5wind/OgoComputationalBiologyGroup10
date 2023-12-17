# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 16:27:55 2023

@author: 20213628
"""
from rdkit import Chem
from rdkit.Contrib.NP_Score import npscorer

def calculate_average_np_score(smiles_list):


    total_np_score = 0.0
    valid_count = 0
    fscore = npscorer.readNPModel()

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Calculate NP score using the sascorer function
            np_score, confidence = npscorer.scoreMolWConfidence(mol, fscore)
            total_np_score += np_score
            valid_count += 1


    if valid_count == 0:
        return 0.0

    average_np_score = total_np_score / valid_count
    return average_np_score

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "CCO", "C1=CC=CC=C1"]

average_np_score = calculate_average_np_score(smiles_list_to_check)

print(average_np_score)



