# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 10:59:50 2023

@author: 20192032
"""

import re
from Elements import elements, lc_elements
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def two_letter_elements(smiles_dataset):
    two_letter_tokens = set()
    pattern_two = re.compile(r'[A-Z][a-z]|[a-z]{2}') #capital+small & small+small
    
    for smiles in smiles_dataset:
        matches_two = pattern_two.findall(smiles)

        # Collect chemically valid two-letter tokens
        two_letter_tokens.update(token for token in matches_two if token in elements)

    #Replace chemically valid two-letter tokens with 'Z' in the dataset
    for i, smiles in enumerate(smiles_dataset):
        for token in two_letter_tokens:
            smiles_dataset[i] = smiles_dataset[i].replace(token, 'Z')
    
    return smiles_dataset, two_letter_tokens


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
    pattern_one = re.compile('[A-Za-z0-9%\-+=#\(\)\[\]]') #hoofdletters, kleine letters + symbolen SMILES
    
    
    for smiles in smiles_dataset:
        matches_one = pattern_one.findall(smiles)

        one_letter_tokens.update(matches_one)
    
    one_letter_tokens = one_letter_tokens.intersection(lc_elements) #removes all tokens not found in lc_elements
    
    char_list = list(one_letter_tokens)
    
    char_to_int = {char: index for index, char in enumerate(char_list)}
    int_to_char = {index: char for index, char in enumerate(char_list)}
    
    char_to_int.update({"E" : len(char_to_int)}) #adding End character
    char_to_int.update({"G" : len(char_to_int)}) #adding starting character
    char_to_int.update({"Z" : len(char_to_int)})
    
    int_to_char.update({len(int_to_char) : "E"})
    int_to_char.update({len(int_to_char) : "G"})
    int_to_char.update({len(int_to_char) : "Z"})
    
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
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    return {'Hydrogen Bond Donors': hbd, 'Hydrogen Bond Acceptors': hba, 'Molecular Weight': mw, 'LogP': logp}

def visualisation_distributions(data, output_filename):
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
    descriptor_list = []
    for smiles in data['canonical_smiles']:
        descriptor_list.append(lipinski_descriptors(smiles))
    dataframe = pd.DataFrame(descriptor_list)
    data = data.reset_index(drop=True)
    dataframe = dataframe.reset_index(drop=True)
    data_with_descriptors = pd.concat([data, dataframe], axis=1)
    
    descriptor_columns = list(data_with_descriptors.columns)[1:]
    plt.figure(figsize=(15, 8))
    plt.suptitle('Density Distributions of Lipinski Descriptors', fontsize=16, fontweight='bold')
    for i, descriptor in enumerate(descriptor_columns, 1):
        plt.subplot(2, 2, i)
        sns.kdeplot(data_with_descriptors[descriptor], fill=True, color=f'C{i}', lw=2)
        plt.xlabel(descriptor)
        plt.ylabel('Density')
    
    data_with_descriptors.to_csv(output_filename, index=False)

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
    #3D array to store one-hot encoded sequences
    one_hot =  np.zeros((data.shape[0], embed+2, len(encode)),dtype=np.int8) #embedded+2 to account for begin and end characters
    
    #loop through each SMILES sequence in data set and convert.
    for i,smile in enumerate(data):
        one_hot[i,0,encode["G"]] = 1
        for j,c in enumerate(smile):
            one_hot[i,j+1,encode[c]] = 1
        one_hot[i,len(smile)+1,encode["E"]] = 1
        
    #return input+output arrays
    return one_hot[:,0:-1,:], one_hot[:,1:,:]

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
    
    #iterates through each sequence in one-hot encoded array
    for i in range(array_data.shape[0]):
        decoded_seq = ""
        for j in range(array_data.shape[1]):
            char_index = np.argmax(array_data[i, j, :])
            decoded_char = decode[char_index]
            
            #append char to decoded_seq untill it is end char.
            if decoded_char != "E":
                decoded_seq += decoded_char
            else:
                break
    
    #append the fully decoded sequence to the list
    gen_mol.append(decoded_seq)
    
    return gen_mol




