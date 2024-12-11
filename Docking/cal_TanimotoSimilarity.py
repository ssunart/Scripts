###두 분자 집합 간의 분자 유사도 계산(Tanimoto similarity, ECFP4, nBit=2048)
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from tqdm import tqdm  # For progress bar

# Configurable parameters
data_enam_file = 'path_to_enam_file.smi'
data_ref_file = 'path_to_ref_file.csv'
output_dir = './output/'
enam_smiles_col = 'smiles'
ref_smiles_col = 'smiles'
ref_name_col = 'name'
enam_id_col = 'id'

# Load data
data_enam = pd.read_csv(data_enam_file, sep='\t', low_memory=False)
data_ref = pd.read_csv(data_ref_file)

# Iterate through reference molecules
for index, ref_row in tqdm(data_ref.iterrows(), total=data_ref.shape[0], desc="Processing Reference Molecules"):
    ref_smiles = ref_row[ref_smiles_col]
    ref_name = ref_row[ref_name_col]
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    
    if not ref_mol:
        print(f"Invalid SMILES '{ref_smiles}' skipped.")
        continue

    ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2)

    # Store results
    results = []
    for enam_index, enam_row in data_enam.iterrows():
        enam_smiles = enam_row[enam_smiles_col]
        enam_id = enam_row[enam_id_col]
        enam_mol = Chem.MolFromSmiles(enam_smiles)
        
        if enam_mol:
            enam_fp = AllChem.GetMorganFingerprintAsBitVect(enam_mol, 2)
            similarity = DataStructs.TanimotoSimilarity(ref_fp, enam_fp)
            results.append((enam_id, Chem.MolToSmiles(enam_mol), similarity))
        else:
            print(f"Invalid SMILES '{enam_smiles}' skipped.")

    # Save results
    result_df = pd.DataFrame(results, columns=['Enam_CatalogID', 'Enam_SMILES', 'Similarity'])
    result_df = result_df.sort_values(by='Similarity', ascending=False)
    result_df.to_csv(f"{output_dir}/similarity_results_for_{ref_name}.csv", index=False)

    print(f"Saved results for {ref_name}")
