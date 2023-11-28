import pandas as pd
from rdkit import Chem
import warnings

# Read the data
data = pd.read_csv('ChEMBL33.csv')

warnings.filterwarnings('ignore', category=FutureWarning)

# Validate and standardize SMILES
valid_smiles = []
invalid_indices = []
for index, row in data.iterrows():
    mol = Chem.MolFromSmiles(row['canonical_smiles'])
    if mol:
        valid_smiles.append(Chem.MolToSmiles(mol))
    else:
        invalid_indices.append(index)

# Handle invalid data
# You can either drop these rows or investigate further
data.drop(invalid_indices, inplace=True)

# Add standardized SMILES to the DataFrame
data['standardized_smiles'] = valid_smiles

# Save the preprocessed data
data.to_csv('preprocessed_ChEMBL33.csv', index=False)
