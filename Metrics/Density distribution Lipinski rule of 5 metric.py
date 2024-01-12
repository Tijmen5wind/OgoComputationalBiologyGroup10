# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 19:00:27 2023

@author: 20213628
"""
def lipinski_descriptors(smiles):
    """
    Calculate Lipinski descriptors for a molecule represented by its SMILES notation.

    Lipinski descriptors are a set of rules used in drug discovery to evaluate the drug-likeness
    of a molecule. This function takes a SMILES representation of a molecule, computes specific
    descriptors, and returns them as four separate lists.

    Parameters:
    - smiles (str): The SMILES notation of the molecule.

    Returns:
    - hbd Hydrogen Bond Donors
    - hba: Hydrogen Bond Acceptors
    - mw: Molecular Weight
    - logp: Log(P)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None, None  # Return None for invalid SMILES
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    return hbd, hba, mw, logp

def list_distributions_with_canonicalization(smiles_list, output_filename):
    """
    List density distributions of Lipinski descriptors for a list of SMILES strings after canonicalization
    

    Parameters:
    - smiles_list (list): List of SMILES strings.
    Example:
    ```python
    smiles_list_to_check = ["CCO", "C1=CC=CC=C1", "CC(C)(C)C(=O)O", "CCN", "C(CCCCCCCC)CCCCCCCC(=O)O", "Invalid"]
    output_filename = 'subset_50000_descr'
    visualisation_distributions(smiles_list_to_check, output_filename)
    ```
    
    Returns:
    - hbd_list (list): List of Hydrogen Bond Donors for each valid molecule.
    - hba_list (list): List of Hydrogen Bond Acceptors for each valid molecule.
    - mw_list (list): List of Molecular Weights for each valid molecule.
    - logp_list (list): List of LogP values for each valid molecule.
    """
    hbd_list, hba_list, mw_list, logp_list = [], [], [], []
    valid_smiles = []  # Store valid SMILES
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Canonicalize SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            hbd, hba, mw, logp = lipinski_descriptors(canonical_smiles)
            if hbd is not None:
                hbd_list.append(hbd)
                hba_list.append(hba)
                mw_list.append(mw)
                logp_list.append(logp)
                valid_smiles.append(canonical_smiles)

    
    return hbd_list, hba_list, mw_list, logp_list




