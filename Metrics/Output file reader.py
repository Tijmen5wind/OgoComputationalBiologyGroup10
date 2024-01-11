# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 19:29:11 2024

@author: 20213628
"""

import random

def extract_first_words_from_csv(filename):
    """
    Extracts the first word from each row in a CSV file and puts them into a list.

    Parameters:
    - filename (str): Name of the CSV file.

    Returns:
    - first_words_list (list): List of the first words from each row.
    """
    # Read the CSV file
    df = pd.read_csv(filename, header=None)

    # Extract the first word from each row
    first_words_list = df.apply(lambda row: str(row[0]).split()[0], axis=1).tolist()

    return first_words_list

def replace_Z_randomly(smiles):
    """
    Replace 'Z' character in a SMILES string with either 'Cl' or 'Br' randomly.

    Parameters:
    - smiles (str): SMILES string.

    Returns:
    - modified_smiles (str): SMILES string with 'Z' replaced.
    """
    return smiles.replace('Z', random.choice(['Cl', 'Br']))

def replace_Z_characters_in_list_randomly(smiles_list):
    """
    Replace 'Z' characters in each SMILES string within a list with either 'Cl' or 'Br' randomly.

    Parameters:
    - smiles_list (list): List of SMILES strings.

    Returns:
    - modified_smiles_list (list): List of SMILES strings with 'Z' replaced.
    """
    modified_smiles_list = [replace_Z_randomly(smiles) for smiles in smiles_list]
    return modified_smiles_list

def remove_first_character_and_empty_elements_random(phrase_list):
    """
    Removes the first character of each element in the list
    and removes the element if nothing is left.
    Replaces 'Z' randomly with 'Cl' or 'Br'.

    Parameters:
    - phrase_list (list): The list of strings.

    Returns:
    - modified_list (list): The modified list.
    """
    modified_list = []

    for element in phrase_list:
        if len(element) > 1:
            # Add the element to the modified list with the first character removed
            modified_list.append(replace_Z_randomly(element[1:]))
    
    return modified_list

# Example usage
filename = 'generated_smiles_10000.csv' # Replace this with the actual filename

first_words = extract_first_words_from_csv(filename)
modified_list_random = remove_first_character_and_empty_elements_random(first_words)

# Replace 'Z' characters in the SMILES list randomly
modified_smiles_list_random = replace_Z_characters_in_list_randomly(modified_list_random)

# Calculate the validity percentage with the modified SMILES list
percentage_valid_random = validity_percentage(modified_smiles_list_random)

print(f"Validity Percentage (Random): {percentage_valid_random:.2f}%")
