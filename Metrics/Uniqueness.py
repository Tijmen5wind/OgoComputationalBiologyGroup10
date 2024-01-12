# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 15:53:31 2023

@author: 20213628
"""

from rdkit import Chem

def uniqueness_with_canonicalization(smiles_list):
    """
    Count the percentage of unique SMILES in a list after canonicalization.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - uniqueness (float): Percentage of unique SMILES.
    """
    # Create RDKit molecules from SMILES strings
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]

    if not mol_list:
        raise ValueError("No valid molecules in the given SMILES list.")

    # Canonicalize SMILES strings
    canonical_smiles_list = [Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True) for mol in mol_list]

    # Count the percentage of unique SMILES
    unique_smiles_set = set(canonical_smiles_list)
    unique_count = len(unique_smiles_set)
    uniqueness = unique_count / len(canonical_smiles_list)

    return uniqueness

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "CCO", "C1=CC=CC=C1"]

uniqueness = uniqueness_with_canonicalization(smiles_list_to_check)

print(f"Uniqueness: {uniqueness:.2f}")
