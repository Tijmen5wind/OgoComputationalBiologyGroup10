# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 00:25:29 2023

@author: 20213628
"""
from rdkit import Chem

def calculate_novelty_with_canonicalization(smiles_list, reference_dataset):
    """
    Calculate the novelty of SMILES strings based on a reference dataset.

    Parameters:
    - smiles_list (list): List of SMILES strings to check for novelty.
    - reference_dataset (list): Reference dataset of SMILES strings.

    Returns:
    - novelty_percentage (float): Percentage of novel SMILES strings.
    """
    # Create RDKit molecules from SMILES strings
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]

    if not mol_list:
        raise ValueError("No valid molecules in the given SMILES list.")

    # Canonicalize SMILES strings
    canonical_smiles_list = [Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True) for mol in mol_list]

    # Calculate novelty
    novel_smiles = sum(smiles not in reference_dataset for smiles in canonical_smiles_list)
    novelty_percentage = (novel_smiles) / len(smiles_list)

    return novelty_percentage

# Example usage:
reference_dataset = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O"]
smiles_to_check = ["CCO", "C1=CC=CC=C1", "InvalidSMILES", "CC(C)(C)C(=O)O", "CCN", "CCN", "CCO"]

novelty_percentage = calculate_novelty_with_canonicalization(smiles_to_check, reference_dataset)

print(f"Novelty Percentage: {novelty_percentage:.2f}")
