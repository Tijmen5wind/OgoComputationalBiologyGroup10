"""
Implementation of all preprocessing steps
"""

import pandas as pd
import numpy as np
from rdkit import Chem
import os

np.random.seed(1)

class Preprocessor:

    def __init__(self, name):
        # where name is the name of the file 

        # List to store data
        self._data = []

        # If True, check after each function that all duplicates are still removed
        self._duplicates_removed = False

        if os.path.isfile(name + '.csv'):
            self._data = pd.read_csv(name + '.csv', header=None).values[:, 0]
        elif os.path.isfile(name + '.tar.xz'):
            # Skip first line since empty and last line since nan
            self._data = pd.read_csv(name + '.tar.xz', compression='xz', header=None).values[1:-1, 0]

        # Remove empty dimensions
        self._data = np.squeeze(self._data)
        return

    def remove_not_valid(self):
        """Remove all SMILES not accepted by the RDKit using a vectorized approach."""
        smiles_series = pd.Series(self._data)
        valid_smiles = smiles_series.apply(lambda s: Chem.MolFromSmiles(str(s)) is not None)
        self._data = smiles_series[valid_smiles].to_numpy()

    def remove_duplicates(self):
        """Remove all SMILES appearing more than once
        :return:
        """
        self._data = np.unique(self._data)

        # Set flag to always remove duplicated after an operation
        self._duplicates_removed = True
        return

    def remove_stereochem(self):
        """Remove all tokens related to stereochemistry using vectorized string operations."""

        smiles_series = pd.Series(self._data)
        stereochem_tokens = ['/', '@', '\\']
        for token in stereochem_tokens:
            smiles_series = smiles_series.str.replace(token, '')
        self._data = smiles_series.to_numpy()

        # Remove possible created duplicates
        if self._duplicates_removed:
            self.remove_duplicates()


    def remove_token(self, t):
        """Remove token t from all elements of data
        :param t:   token to remove
        :return:
        """
        self._data = np.array([d.replace(t, '') for d in self._data])

        # Remove possible created duplicates
        if self._duplicates_removed:
            self.remove_duplicates()
        return

    def remove_salts(self):
        """Remove all salts by retaining the longest part of the SMILES sequence."""

        smiles_series = pd.Series(self._data)
        # Select longest part of SMILES
        cleaned_smiles = smiles_series.apply(lambda s: max(s.split('.'), key=len))
        self._data = cleaned_smiles.to_numpy()

        # Remove possible deposits
        self.remove_token('.')

        # Remove possible created duplicates
        if self._duplicates_removed:
            self.remove_duplicates()


    def canonicalize(self):
        """Canonicalize all SMILES from data."""
        smiles_series = pd.Series(self._data)
        canonical_smiles = smiles_series.apply(lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(str(s)), isomericSmiles=True, canonical=True))
        self._data = canonical_smiles.to_numpy()
        if self._duplicates_removed:
            self.remove_duplicates()

    def save_data(self, name='data.csv'):
        """Saves data to file"""
        pd.DataFrame(self._data).to_csv(name, header=None, index=None)
        return

    def get_data(self):
        """Returns data"""
        return self._data
