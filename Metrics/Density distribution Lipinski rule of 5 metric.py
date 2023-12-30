# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 19:00:27 2023

@author: 20213628
"""
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import seaborn as sns

def lipinski_descriptors(smiles):
    """
    Calculate Lipinski descriptors for a molecule represented by its SMILES notation.

    Lipinski descriptors are a set of rules used in drug discovery to evaluate the drug-likeness
    of a molecule. This function takes a SMILES representation of a molecule, computes specific
    descriptors, and returns them as a dictionary.

    Parameters:
    - smiles (str): The SMILES notation of the molecule.

    Returns:
    dict: A dictionary containing Lipinski descriptors.
        - 'Hydrogen Bond Donors': Number of hydrogen bond donors.
        - 'Hydrogen Bond Acceptors': Number of hydrogen bond acceptors.
        - 'Molecular Weight': Molecular weight of the molecule.
        - 'LogP': The logarithm of the partition coefficient between octanol and water.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None  # Return None for invalid SMILES
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    return {'Hydrogen Bond Donors': hbd, 'Hydrogen Bond Acceptors': hba, 'Molecular Weight': mw, 'LogP': logp}

def visualisation_distributions_with_canonicalization(smiles_list, output_filename):
    """
    Visualize density distributions of Lipinski descriptors for a list of SMILES strings after canonicalization
    and save the results to a CSV file.

    This function calculates Lipinski descriptors for each valid molecule in the provided list of SMILES strings,
    canonicalizes the SMILES strings, creates a DataFrame with the descriptors, and visualizes the density distributions
    of the descriptors. The resulting plots are saved to a file, and the DataFrame with descriptors is saved to a CSV file.

    Parameters:
    - smiles_list (list): List of SMILES strings.
    - output_filename (str): The filename for saving the visualizations and the DataFrame with descriptors.

    Example:
    ```python
    smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "C(CCCCCCCC)CCCCCCCC(=O)O", "Invalid"]
    output_filename = 'subset_50000_descr.csv'
    visualisation_distributions(smiles_list_to_check, output_filename)
    ```
    Output:
    - Saved visualizations in the current directory.
    - Saved DataFrame with descriptors as 'subset_50000_descr.csv'.
    """
    descriptor_list = []
    valid_smiles = []  # Store valid SMILES
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Canonicalize SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            descriptor = lipinski_descriptors(canonical_smiles)
            if descriptor is not None:
                descriptor_list.append(descriptor)
                valid_smiles.append(canonical_smiles)
    
    dataframe = pd.DataFrame(descriptor_list)
    
    descriptor_columns = list(dataframe.columns)
    plt.figure(figsize=(15, 8))
    plt.suptitle('Density Distributions of Lipinski Descriptors', fontsize=16, fontweight='bold')
    for i, descriptor in enumerate(descriptor_columns, 1):
        plt.subplot(2, 2, i)
        sns.kdeplot(dataframe[descriptor], fill=True, color=f'C{i}', lw=2)
        plt.xlabel(descriptor)
        plt.ylabel('Density')
    
    dataframe.to_csv(output_filename, index=False)
    return valid_smiles

# Example usage:
smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "C(CCCCCCCC)CCCCCCCC(=O)O", "Invalid"]
output_filename = 'subset_50000_descr.csv'
valid_canonical_smiles = visualisation_distributions_with_canonicalization(smiles_list_to_check, output_filename)




