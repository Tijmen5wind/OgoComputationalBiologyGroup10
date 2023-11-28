import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Load the preprocessed data
data = pd.read_csv('preprocessed_ChEMBL33.csv')

# Function to calculate Lipinski descriptors
def lipinski_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    return hbd, hba, mw, logp

# Apply the function to each row in the DataFrame
data[['HBD', 'HBA', 'MW', 'LogP']] = data.apply(
    lambda row: lipinski_descriptors(row['standardized_smiles']), axis=1, result_type='expand'
)

# Save the data with descriptors
data.to_csv('ChEMBL33_with_descriptors.csv', index=False)
