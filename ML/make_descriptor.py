###도킹 데이터 처리 & molecular descriptor를 계산 및 저장

#!/bin/bash
#SBATCH --partition=rome_short    # cpu1
# #SBATCH --partition=skylake_short # cpu[2-7]
#SBATCH --ntasks=6
#SBATCH --mem=0
#SBATCH --job-name=sort_adGPU_p1
#SBATCH --output=sort_adGPU01.out

module load compiler/gcc-9.4.0 cuda/11.3

!pip install rdkit

#Module load for dataprep#

import pandas as pd
import numpy as np
from tqdm import tqdm

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sklearn
from tqdm import tqdm # to show a progress bar.

#Load docking data#


df_raw = pd.read_csv("refined_DL_gpuR.csv") #csv file load


###일부만 가져오기###

#df_raw.describe()

#df_sample = df_raw.sample(frac=0.01) #랜덤하게 10% 가져오기

#####################

#Data 전처리: canonical smiles#

# SMILES, 원하는 score 추출하여 canonical smiles 전환 후 새 csv file(dock_canonical.csv) 셋업
fout = open("dock_canonical.csv", "w")
fout.write("CatalogID,SMILES,AD_GPU\n")

idx = 0
for _, row in df_sample.iterrows():
    CID = row['CatalogID']
    smiles = row['SMILES']
    score = row['AD_GPU']
    if pd.isna(score):  # if  score is NaN, continue to next row
        continue
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:  # if not a valid SMILES.
        continue
    new_smi = Chem.MolToSmiles(mol, canonical=True)

    fout.write(f"{CID},{new_smi},{score}\n")
    idx += 1

print(f"Total number of converted molecules: {idx}") #정상적으로 전환된 분자 수 출력
fout.close()

data = pd.read_csv("dock_canonical.csv")


"""data type check"""

# 숫자가 아닌 값이 있는지 확인합니다.
is_non_numeric = pd.to_numeric(data['AD_GPU'], errors='coerce').isna()

# 숫자가 아닌 값이 있는 행을 삭제합니다.
df_a = data[~is_non_numeric]

# 인덱스를 재설정합니다. drop=True로 설정하면 이전 인덱스를 삭제합니다.
df_a = df_a.reset_index(drop=True)


# 변경 사항을 새 csv 파일에 저장합니다.
df_a.to_csv('cleaned.csv', index=False)

df_a['AD_GPU'] = df_a['AD_GPU'].astype(float)

df = pd.DataFrame(df_a)

score_df=df["AD_GPU"]

##Molecular descriptor: MACCS key 추가##

from rdkit import Chem
from rdkit.Chem import MACCSkeys

# 분자들의 리스트를 생성
mols = [Chem.MolFromSmiles(x) for x in df["SMILES"]]

# MACCS 키 생성하는 함수를 정의
def generate_MACCS(mol):
    maccs_key = MACCSkeys.GenMACCSKeys(mol)
    return list(map(int, maccs_key.ToBitString()))  # binary string to int list

# 각 분자에 대해 MACCS 키를 생성하고, 이를 리스트에 추가
maccs_list = [generate_MACCS(m) for m in mols]

# 리스트를 데이터프레임으로 변환
maccs_df = pd.DataFrame(maccs_list)


##Molecular descriptor: ECFP 추가##

mols = [Chem.MolFromSmiles(x) for x in df["SMILES"]]

radius = 2; nbits = 1024 # radius=1 ECFP2 / radius=2 ECFP4
fp_list = []
for m in tqdm(mols):
  fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits) # if radius = 2, then fp is an ECFP4
  fp_list.append(fp.ToList())

fp_df = pd.DataFrame(fp_list)

###Mordred descriptor: Mordred 추가##

!pip install mordred

from rdkit import Chem
from mordred import Calculator, descriptors

# create descriptor calculator with all descriptors
calc = Calculator(descriptors, ignore_3D=True)


mordred_df = calc.pandas(mols)

object_dtype = [ idx for idx, i in enumerate(mordred_df.dtypes) if i == 'object' ]
new_mordred_df = mordred_df.drop(mordred_df.columns[object_dtype], axis=1)


###descriptor 합치기###

###dataframe 합치기
merged_df = pd.concat([df_a['SMILES'], score_df, fp_df, maccs_df], axis=1)  ##new_mordred_df
merged_df.to_csv('mergedDL_desc.csv', index=False)
merged_df.to_pickle('mergedDL_desc.pkl')

"""1 단계: 데이터준비 완료 - 분기점

"""

dock_data = pd.read_pickle('merged5K_desc.pkl')

dock_data

"""특정 열 숫자확인"""

#cols_to_check = merged_df.columns.difference(['CID', 'smiles'])

#is_numeric = merged_df[cols_to_check].applymap(np.isreal).all()
##모든 열 확
is_numeric_all = dock_data.applymap(np.isreal).all()
print(is_numeric_all)
