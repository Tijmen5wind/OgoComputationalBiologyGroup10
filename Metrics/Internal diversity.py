# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 19:07:48 2023

@author: 20213628
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np

def calculate_ecfp(mol_list, radius=3, length=1024):
    """
    Calculate ECFP fingerprints for a list of molecules.

    Parameters:
    - mol_list (list): List of RDKit molecule objects.
    - radius (int): ECFP fingerprint radius.
    - length (int): ECFP fingerprint length.

    Returns:
    - fingerprints (list): List of ECFP fingerprints.
    """
    fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length) for mol in mol_list]
    return fingerprints

def calculate_tanimoto_similarity(fp1, fp2):
    """
    Calculate Tanimoto similarity between two fingerprints.

    Parameters:
    - fp1, fp2 (rdkit.DataStructs.cDataStructs.ExplicitBitVect): ECFP fingerprints.

    Returns:
    - similarity (float): Tanimoto similarity.
    """
    return TanimotoSimilarity(fp1, fp2)

def calculate_internal_diversity_with_canonicalization(smiles_list):
    """
    Calculate internal diversity for a list of SMILES strings after canonicalization.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - internal_diversity (float): Internal diversity.
    """
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]

    if not mol_list:
        return None

    # Canonicalize SMILES
    canonical_smiles_list = [Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True) for mol in mol_list]

    fingerprints = calculate_ecfp(mol_list)

    similarities = []
    for i in range(len(canonical_smiles_list)):
        for j in range(i + 1, len(canonical_smiles_list)):
            similarity = calculate_tanimoto_similarity(fingerprints[i], fingerprints[j])
            similarities.append(similarity)

    internal_diversity = np.mean(similarities)

    return internal_diversity

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "C(CCCCCCCC)CCCCCCCC(=O)O"]

internal_diversity = calculate_internal_diversity_with_canonicalization(smiles_list_to_check)

if internal_diversity is not None:
    print('Internal diversity:', internal_diversity)
else:
    print("No valid molecules in the given SMILES list.")

