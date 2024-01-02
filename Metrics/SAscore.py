# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 16:13:35 2023

@author: 20213628
"""

from rdkit import Chem
from rdkit.Contrib.SA_Score import sascorer
import seaborn as sns
import matplotlib.pyplot as plt

def plot_sa_scores_distribution_with_canonicalization(smiles_list):
    """
    Plot the distribution of Synthetic Accessibility (SA) scores for a list of molecules after canonicalization.

    Parameters:
    - smiles_list (list): List of SMILES strings.
    """

    sa_scores = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Canonicalize SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            # Calculate SA score using the sascorer function
            sa_score = sascorer.calculateScore(Chem.MolFromSmiles(canonical_smiles))
            sa_scores.append(sa_score)

    # Plot the distribution using Seaborn
    sns.kdeplot(sa_scores, color='blue', fill=True)
    plt.title('Distribution of SA Scores')
    plt.xlabel('SA Score')
    plt.ylabel('Density')
    plt.show()

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "CCO", "C1=CC=CC=C1", "HOI"]
plot_sa_scores_distribution_with_canonicalization(smiles_list_to_check)


