###CSV 파일에서 SMILES 데이터(화학 구조를 문자열로 표현)와 해당 ID를 추출

import pandas as pd

# Load the CSV file
file_path = 'merged_selections_CACHE4.csv'  # Replace with your file path
data = pd.read_csv(file_path)

# Extracting only the 'Smiles' column
smiles_only = data[['Smiles','Merged_ID']]

# Save this as a new text file with each SMILES string on a new line
smiles_file_path = '1800ligands.smi'  # Define your output file path
smiles_only.to_csv(smiles_file_path, index=False, header=False)
