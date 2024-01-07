# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 10:59:50 2023

@author: 20192032
"""

import re
from Elements import lc_elements
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from Elements import elements
import matplotlib.pyplot as plt
import seaborn as sns


def data_filtering(smiles_dataset):
    """
    Filters and modifies a dataset of SMILES (Simplified Molecular Input Line Entry System) notations.

    Args:
    - smiles_dataset (list): A list of SMILES strings.

    Returns:
    - modified_dataset (list): A list of modified SMILES strings where 'Br' and 'Cl' are replaced with 'Z'.
    - two_letter_tokens (set): A set of unique two-letter tokens found in the dataset.
    """

    # Initialize a set to store unique two-letter tokens
    two_letter_tokens = set()

    # Compile a regex pattern to find two-letter tokens
    pattern_two = re.compile(r'\b[A-Z][a-z]?\b|\b[a-z]{2}\b')

    # Iterate over the SMILES dataset
    for smiles in smiles_dataset:
        # Find all matches for the pattern in each SMILES string
        matches_two = pattern_two.findall(smiles)

        # Update the set with valid two-letter tokens found in SMILES
        two_letter_tokens.update(token for token in matches_two if token in elements)

    # Remove specific elements ('Br' and 'Cl') from the set of tokens
    two_letter_tokens.discard('Br')
    two_letter_tokens.discard('Cl')

    # Filter out SMILES strings containing any of the two-letter tokens
    filtered_smiles = [
        smiles for smiles in smiles_dataset
        if not any(element in smiles for element in two_letter_tokens)
    ]

    # Define tokens to be replaced
    replace = {'Br', 'Cl'}
    modified_dataset = []

    # Modify the filtered dataset by replacing specified tokens with 'Z'
    for smiles in filtered_smiles:
        modified_smiles = smiles
        for token in replace:
            modified_smiles = modified_smiles.replace(token, 'Z')
        modified_dataset.append(modified_smiles)

    return modified_dataset, two_letter_tokens

def create_dictionary(smiles_dataset):
    """
    Create dictionaries for tokenizing SMILES strings.

    This function takes a dataset of SMILES strings and creates dictionaries
    for tokenizing the strings into integer sequences. It extracts symbols, one-letter
    and two-letter chemical element tokens from the SMILES strings, checks if they are chemically valid
    and assigns unique integer indices to each token.

    Parameters:
    - smiles_dataset (list): A list of SMILES strings.

    Returns:
    tuple: A tuple of dictionaries.
        - char_to_int (dict): A dictionary mapping characters to integer indices.
        - int_to_char (dict): A dictionary mapping integer indices to characters.
    
    imports from Elements module:
    - `elements`: list of chemical valid two-letter chemical elements.
    - `lc_elements`: list of one-letter chemical elements (both lower and uppercase) and valid SMILES symbols. 
    """
    one_letter_tokens = set()
    pattern_one = re.compile('Br|Cl|[A-Za-z0-9%\-+=#\(\)\[\]]')  # hoofdletters, kleine letters + symbolen SMILES

    for smiles in smiles_dataset:
        matches_one = pattern_one.findall(smiles)

        one_letter_tokens.update(matches_one)

    one_letter_tokens = one_letter_tokens.intersection(lc_elements)  # removes all tokens not found in lc_elements

    char_list = list(one_letter_tokens)

    char_to_int = {char: index for index, char in enumerate(char_list)}
    int_to_char = {index: char for index, char in enumerate(char_list)}

    char_to_int.update({"E": len(char_to_int)})  # adding End character
    char_to_int.update({"G": len(char_to_int)})  # adding starting character
    char_to_int.update({"Z": len(char_to_int)})

    int_to_char.update({len(int_to_char): "E"})
    int_to_char.update({len(int_to_char): "G"})
    int_to_char.update({len(int_to_char): "Z"})

    return char_to_int, int_to_char


def subset_splitting(datafile, num_subsets=6, step_size=25000):
    """
    Generate subset sizes for splitting a dataset into multiple subsets.

    This function calculates a list of subset sizes for splitting a dataset
    into a specified number of subsets with increasing sample sizes. The sample
    sizes are determined based on an initial size (50000) plus multiples of a
    specified step size. The function excludes subset sizes that exceed the
    total size of the dataset.

    Parameters:
    - datafile (list or array-like): The dataset to be split.
    - num_subsets (int, optional): The number of subsets to generate. Default is 6.
    - step_size (int, optional): The step size to increment the sample sizes. Default is 25000.

    Returns:
    list: A list of subset sizes.
    """
    sample_sizes = [50000 + (i * step_size) for i in range(num_subsets)]
    sample_sizes = [size for size in sample_sizes if size < len(datafile)]
    return sample_sizes


def subset_creation(smiles_list, sample_sizes):
    """
    Create subsets of a dataset with specified sample sizes.

    This function generates subsets of a given dataset by randomly sampling
    data points with replacement. The sample sizes are provided as a list,
    and subsets are created for each specified group size. The function returns
    a dictionary where the keys are the group sizes, and the values are the
    corresponding subsets.

    Parameters:
    - datafile (DataFrame or array-like): The dataset from which to create subsets.
    - sample_sizes (list): A list of integers specifying the sample sizes for each subset.

    Returns:
    dict: A dictionary where keys are group sizes, and values are corresponding subsets.
    """
    df = pd.DataFrame({'canonical_smiles': smiles_list})
    subsets = {}  # Dictionary to store subsets with group size as key
    for groupsize in sample_sizes:
        subgroup = df.sample(groupsize, replace=True, random_state=42)
        subsets[groupsize] = subgroup
    return subsets


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
    if mol is not None:

        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        return {'HBD': hbd, 'HBA': hba, 'MW': mw, 'LogP': logp}

    else:
        # Return a default dictionary with NaN values
        return {'HBD': float('nan'), 'HBA': float('nan'), 'MW': float('nan'), 'LogP': float('nan')}


def visualisation_distributions(subset_dict, output_filename):
    """
    Visualize density distributions of Lipinski descriptors and save the results to a CSV file.

    This function calculates Lipinski descriptors for each molecule in the provided dataset,
    creates a DataFrame with the descriptors, and visualizes the density distributions of the descriptors.
    The resulting plots are saved to a file, and the DataFrame with descriptors is saved to a CSV file.

    Parameters:
    - data (DataFrame): The input dataset containing the 'canonical_smiles' column.
    - output_filename (str): The filename for saving the visualizations and the DataFrame with descriptors.

    Returns:
    None

    Example:
    ```python
    data_example = pd.DataFrame({'canonical_smiles': ['CCO', 'CCC', 'CCCC'],
                                 'other_column': [1, 2, 3]})
    visualisation_distributions(data_example, 'lipinski_distributions.csv')
    ```
    Output:
    - Saved visualizations in the current directory.
    - Saved DataFrame with descriptors as 'lipinski_distributions.csv'.
    """
    all_descriptors = []

    # Collecting descriptors for all subsets
    for key, data in subset_dict.items():
        descriptor_list = [lipinski_descriptors(smiles) for smiles in data['canonical_smiles']]
        for descriptor in descriptor_list:
            descriptor['Subset Size'] = key  # Add subset size for identification
        all_descriptors.extend(descriptor_list)

    # Creating a combined DataFrame
    combined_data = pd.DataFrame(all_descriptors)

    # Set the figure size and plot the distributions
    plt.figure(figsize=(15, 10))

    for i, descriptor in enumerate(['HBD', 'HBA', 'MW', 'LogP'], 1):
        plt.subplot(2, 2, i)
        for subset_size in subset_dict.keys():
            subset_data = combined_data[combined_data['Subset Size'] == subset_size]
            sns.kdeplot(subset_data[descriptor], label=f'Subset {subset_size}', fill=True, lw=2)
        plt.xlabel(descriptor)
        plt.ylabel('Density')
        plt.title(f'Density Distribution of {descriptor}')
        plt.legend()

    plt.suptitle('Density Distributions of Lipinski Descriptors Across Subsets', fontsize=16, fontweight='bold')
    plt.savefig(output_filename)

def encoder(data, encode):
    """
    Encodes SMILES sequences into one-hot encoding.

    Parameters:
        data (numpy.ndarray): Array containing SMILES sequences.
        char_to_int (dict): Mapping from characters to integer indices.
        max_seq (int): Length of the longest sequence in the dataset for padding.

    Returns:
        tuple: Tuple containing input and output one-hot encoded arrays.
    """
    embed = max([len(seq) for seq in data])
    # 3D array to store one-hot encoded sequences
    one_hot = np.zeros((data.shape[0], embed + 2, len(encode)),
                       dtype=np.int8)  # embedded+2 to account for begin and end characters

    # loop through each SMILES sequence in data set and convert.
    for i, smile in enumerate(data):
        one_hot[i, 0, encode["G"]] = 1
        for j, c in enumerate(smile):
            one_hot[i, j + 1, encode[c]] = 1
        one_hot[i, len(smile) + 1, encode["E"]] = 1

    # return input+output arrays
    return one_hot[:, 0:-1, :], one_hot[:, 1:, :]


def decoder(array_data, decode):
    """
    Decodes one-hot encoded sequences into SMILES strings.

    Parameters:
        array_data (numpy.ndarray): 3D array with 2D one-hot encoded array to be decoded.
        int_to_char (dict): Mapping from integer indices to characters.

    Returns:
        list: list of decoded SMILES strings.
    """
    gen_mol = []

    # iterates through each sequence in one-hot encoded array
    for i in range(array_data.shape[0]):
        decoded_seq = ""
        for j in range(array_data.shape[1]):
            char_index = np.argmax(array_data[i, j, :])
            decoded_char = decode[char_index]

            # append char to decoded_seq untill it is end char.
            if decoded_char != "E":
                decoded_seq += decoded_char
            else:
                break

    # append the fully decoded sequence to the list
    gen_mol.append(decoded_seq)

    return gen_mol
