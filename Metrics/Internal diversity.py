# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 19:07:48 2023

@author: 20213628
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
import itertools
import random

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

def calculate_internal_diversity(smiles_list, sample_size=10000):
    """
    Calculate internal diversity for a list of SMILES strings.

    Parameters:
    - smiles_list (list): List of SMILES strings.
    - sample_size (int): Number of random pairs to analyze.

    Returns:
    - internal_diversity (float): Internal diversity.
    """
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]

    if not mol_list:
        return None

    fingerprints = calculate_ecfp(mol_list)

    # Generate a random sample of pairs
    pairs_to_compare = random.sample(list(itertools.combinations(range(len(mol_list)), 2)), min(sample_size, len(mol_list)*(len(mol_list)-1)//2))

    similarities = []
    for pair in pairs_to_compare:
        i, j = pair
        similarity = calculate_tanimoto_similarity(fingerprints[i], fingerprints[j])
        similarities.append(similarity)

    internal_diversity = np.mean(similarities)

    return internal_diversity

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "C(CCCCCCCC)CCCCCCCC(=O)O"]

internal_diversity = calculate_internal_diversity(smiles_list_to_check)

if internal_diversity is not None:
    print('internal diversity:',internal_diversity)
else:
    print("No valid molecules in the given SMILES list.")
