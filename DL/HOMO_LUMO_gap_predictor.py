#Load rdkit & version check
import rdkit
rdkit.__version__
import rdkit.Chem

###Reading the data### 
#reference: https://pubs.acs.org/doi/10.1021/acs.jcim.7b00083
#PubChemQC database contains about 200 million TD-DFT calculation results using the GAMESS QM package.
#The DB contains the following information
#Optimized molecular structures
#HOMO-LUMO gap
#Excitation state information, such as oscillator strength and energies of excited states
#osciilator strength ~ absoption coefficient

import pandas as pd
pubchem = pd.read_csv("https://www.dropbox.com/s/s9uhxw06z42gs8b/PubchemQC_subset_HOMO-LUMO_and_OS.csv?dl=1", on_bad_lines='skip')

#10만 개 data 사용
from tqdm import tqdm # progressive bar 표시를 위해서... 
gap_list = []
os_list = []
smi_list = []
ii = 0
max_mol = 100000 # 일단 10만개의 데이터를 사용하자... 
for gap, os, smi in zip(tqdm(pubchem["HOMO-LUMO_gap(eV)"]), pubchem["Oscillator_Strength"], pubchem["SMILES"]):
  gap = float(gap)
  os = float(os)
  gap_list.append(gap) # HOMO-LUMO gap
  os_list.append(os)  # Oscillator strength 값
  smi_list.append(smi) # 분자의 SMILES 표현식
  ii += 1
  if ii >= max_mol: 
    break
  else:
    continue

print(len(gap_list))
print(gap_list)
print(os_list)
print(smi_list[:10])

###Converting a molecule into an extended circular fingerprint (ECFP)
from rdkit.Chem import AllChem
mol = rdkit.Chem.MolFromSmiles(smi_list[0]) # SMILES -> Mol-type 변수로 읽어들인다. 
vec = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 1024) # ECFP2를 이용. 1024개의 0/1로 표현하겠다.

smi_list[0]
print(vec)
print(vec.ToList())

#RDKit을 이용해서 10만개 분자를 ECFP로 변환
X_ECFP = []
for smi in tqdm(smi_list):
  mol = rdkit.Chem.MolFromSmiles(smi)
  vec = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 1024)
  X_ECFP.append(vec.ToList()) # rdkit 최근 버젼.
  #X_ECFP.append([int(x) for x in vec.ToBitString()]) # rdkit 구 버젼 (2021 이전 버젼...)

print(len(X_ECFP))
print(X_ECFP[:2])
Y = gap_list # 예측하고자하는 목표 물성: HOMO-LUMO gap...

###Spliting training and test set###
#The ratio of training, validation, and test (holdout) data is a little bit arbitrary.
#The most common ratios are 6:2:2, 8:1:1, and 7:1:2 etc.

#Scikit-learn 로드
import sklearn
from sklearn.model_selection import train_test_split


