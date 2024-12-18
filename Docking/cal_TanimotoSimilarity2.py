### RDKit 라이브러리를 사용하며, 각 ref_setA의 분자를 기준으로 모든 REAL_scr 분자와의 타니모토 유사도를 계산
#Module
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs


data_enam = pd.read_csv('/home/soowon/project/07_PPAR/REAL_scr/canonical/HAC22_23/smi_14_re.smi',sep='\t',low_memory=False)
data_ref =  pd.read_csv('/home/soowon/project/07_PPAR/ref_setA.csv')


# data_ref의 모든 SMILES에 대해 반복
for index, ref_row in data_ref.iterrows():
    ref_smiles = ref_row['smiles']
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    if ref_mol:
        ref_mol = Chem.MolToSmiles(ref_mol, canonical=True)  # Canonical SMILES로 변환
        ref_mol = Chem.MolFromSmiles(ref_mol)  # Canonical SMILES로부터 분자 객체 생성
        ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2)
    else:
        print(f"Invalid SMILES '{ref_smiles}' skipped.")
        continue

    # 결과를 저장할 리스트
    results = []

    # data_enam의 모든 SMILES에 대해 반복
    for enam_index, enam_row in data_enam.iterrows():
        enam_smiles = enam_row[0]  # 열 이름 확인 필요
        enam_mol = Chem.MolFromSmiles(enam_smiles)
        if enam_mol:
            enam_mol = Chem.MolToSmiles(enam_mol, canonical=True)  # Canonical SMILES로 변환
            enam_mol = Chem.MolFromSmiles(enam_mol)  # Canonical SMILES로부터 분자 객체 생성
            enam_fp = AllChem.GetMorganFingerprintAsBitVect(enam_mol, 2)
            similarity = DataStructs.TanimotoSimilarity(ref_fp, enam_fp)

            # 결과 리스트에 추가
            results.append((enam_row[1], Chem.MolToSmiles(enam_mol), similarity))  # Canonical SMILES 사용
        else:
            print(f"Invalid SMILES '{enam_smiles}' skipped.")

    # 결과를 DataFrame으로 변환하고 유사도 순으로 정렬
    result_df = pd.DataFrame(results, columns=['Enam_CatalogID', 'Enam_SMILES', 'Similarity'])
    result_df = result_df.sort_values(by='Similarity', ascending=False)

    # 각 리간드별로 결과 파일 저장
    file_name = f"/home/soowon/project/07_PPAR/REAL_scr/canonical/HAC22_23/sim14R_{ref_row['name']}.csv"
    result_df.to_csv(file_name, index=False)

    print(f"Saved results for {ref_row['name']} to {file_name}")
