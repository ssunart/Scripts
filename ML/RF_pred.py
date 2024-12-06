#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


df_raw = pd.read_csv("refined_DL_gpuR.csv") #csv file load
df_raw  ##Check the file

df_raw.describe()


# In[3]:


df_sample = df_raw.sample(frac=0.01) #랜덤하게 10% 가져오기

"""startified sampling 시도"""

# 'AD_GPU' 값을 10개의 구간으로 나눕니다.
df_raw['AD_GPU_bin'] = pd.cut(df_raw['AD_GPU'], bins=10, labels=False)

# Stratified sampling을 수행합니다.
from sklearn.model_selection import train_test_split

_, df_sample = train_test_split(df_raw, test_size=0.1, stratify=df_raw['AD_GPU_bin'], random_state=42)

# 'AD_GPU_bin' 컬럼을 삭제합니다.
df_sample = df_sample.drop(columns='AD_GPU_bin')

df_sample


# In[4]:


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
data ##Data check


# In[5]:


# 숫자가 아닌 값이 있는지 확인합니다.
is_non_numeric = pd.to_numeric(data['AD_GPU'], errors='coerce').isna()

# 숫자가 아닌 값이 있는 행을 삭제합니다.
df_a = data[~is_non_numeric]

# 인덱스를 재설정합니다. drop=True로 설정하면 이전 인덱스를 삭제합니다.
df_a = df_a.reset_index(drop=True)

# 결과를 확인합니다.
print(df_a)

# 변경 사항을 새 csv 파일에 저장합니다.
df_a.to_csv('cleaned.csv', index=False)

df_a['AD_GPU'] = df_a['AD_GPU'].astype(float)
df_a['AD_GPU'].describe()
df_a.hist()

df = pd.DataFrame(df_a)


# In[6]:


score_df=df["AD_GPU"]


# In[ ]:



"""MACCS key 추가"""

from rdkit import Chem
from rdkit.Chem import MACCSkeys

# 분자들의 리스트를 생성합니다.
mols = [Chem.MolFromSmiles(x) for x in df["SMILES"]]

# MACCS 키를 생성하는 함수를 정의합니다.
def generate_MACCS(mol):
    maccs_key = MACCSkeys.GenMACCSKeys(mol)
    return list(map(int, maccs_key.ToBitString()))  # binary string to int list

# 각 분자에 대해 MACCS 키를 생성하고, 이를 리스트에 추가합니다.
maccs_list = [generate_MACCS(m) for m in mols]

# 리스트를 데이터프레임으로 변환하고 확인합니다.
maccs_df = pd.DataFrame(maccs_list)
maccs_df  # 데이터 확인


# In[15]:


mols = [Chem.MolFromSmiles(x) for x in df["SMILES"]]
radius = 1; nbits = 2048 # radius=1 ECFP2 / radius=2 ECFP4
fp_list = []
for m in tqdm(mols):
  fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits) # if radius = 2, then fp is an ECFP4
  fp_list.append(fp.ToList())


# In[16]:


fp_df = pd.DataFrame(fp_list)
fp_df


# In[ ]:


import numpy as np
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.GraphDescriptors import (BalabanJ, BertzCT, Chi0, Chi0n, Chi0v, Chi1,
                                         Chi1n, Chi1v, Chi2n, Chi2v, Chi3n, Chi3v, Chi4n, Chi4v,
                                         HallKierAlpha, Ipc, Kappa1, Kappa2, Kappa3)

from rdkit.Chem.EState.EState_VSA import (EState_VSA1, EState_VSA10, EState_VSA11, EState_VSA2, EState_VSA3,
                                          EState_VSA4, EState_VSA5, EState_VSA6, EState_VSA7, EState_VSA8, EState_VSA9,
                                          VSA_EState1, VSA_EState10, VSA_EState2, VSA_EState3, VSA_EState4, VSA_EState5,
                                          VSA_EState6, VSA_EState7, VSA_EState8, VSA_EState9,)

from rdkit.Chem.Descriptors import (ExactMolWt, MolWt, HeavyAtomMolWt, MaxAbsPartialCharge, MinPartialCharge,
                                    MaxPartialCharge, MinAbsPartialCharge, NumRadicalElectrons, NumValenceElectrons)
from rdkit.Chem.EState.EState import (MaxAbsEStateIndex, MaxEStateIndex, MinAbsEStateIndex, MinEStateIndex,)
from rdkit.Chem.Lipinski import (FractionCSP3, HeavyAtomCount, NHOHCount, NOCount, NumAliphaticCarbocycles,
                                 NumAliphaticHeterocycles, NumAliphaticRings, NumAromaticCarbocycles, NumAromaticHeterocycles,
                                 NumAromaticRings, NumHAcceptors, NumHDonors, NumHeteroatoms, RingCount,
                                 NumRotatableBonds, NumSaturatedCarbocycles, NumSaturatedHeterocycles, NumSaturatedRings,)
from rdkit.Chem.Crippen import (MolLogP, MolMR,)
from rdkit.Chem.MolSurf import (LabuteASA, PEOE_VSA1, PEOE_VSA10, PEOE_VSA11, PEOE_VSA12, PEOE_VSA13, PEOE_VSA14,
                                PEOE_VSA2, PEOE_VSA3,PEOE_VSA4, PEOE_VSA5, PEOE_VSA6, PEOE_VSA7, PEOE_VSA8, PEOE_VSA9,
                                SMR_VSA1, SMR_VSA10, SMR_VSA2, SMR_VSA3, SMR_VSA4, SMR_VSA5, SMR_VSA6,
                                SMR_VSA7, SMR_VSA8, SMR_VSA9, SlogP_VSA1, SlogP_VSA10, SlogP_VSA11, SlogP_VSA12,
                                SlogP_VSA2, SlogP_VSA3,SlogP_VSA4, SlogP_VSA5, SlogP_VSA6, SlogP_VSA7, SlogP_VSA8,
                                SlogP_VSA9, TPSA,)
from rdkit.Chem.Fragments import (fr_Al_COO, fr_Al_OH, fr_Al_OH_noTert, fr_ArN, fr_Ar_COO, fr_Ar_N, fr_Ar_NH,
 fr_Ar_OH, fr_COO, fr_COO2, fr_C_O, fr_C_O_noCOO, fr_C_S, fr_HOCCN, fr_Imine, fr_NH0, fr_NH1,
 fr_NH2, fr_N_O, fr_Ndealkylation1, fr_Ndealkylation2, fr_Nhpyrrole, fr_SH, fr_aldehyde, fr_alkyl_carbamate,
 fr_alkyl_halide, fr_allylic_oxid, fr_amide, fr_amidine, fr_aniline, fr_aryl_methyl, fr_azide, fr_azo, fr_barbitur,
 fr_benzene, fr_benzodiazepine, fr_bicyclic, fr_diazo, fr_dihydropyridine, fr_epoxide, fr_ester, fr_ether, fr_furan,
 fr_guanido, fr_halogen, fr_hdrzine, fr_hdrzone, fr_imidazole, fr_imide, fr_isocyan, fr_isothiocyan, fr_ketone,
 fr_ketone_Topliss, fr_lactam, fr_lactone, fr_methoxy, fr_morpholine, fr_nitrile, fr_nitro, fr_nitro_arom,
 fr_nitro_arom_nonortho, fr_nitroso, fr_oxazole, fr_oxime, fr_para_hydroxylation, fr_phenol,
 fr_phenol_noOrthoHbond, fr_phos_acid, fr_phos_ester, fr_piperdine, fr_piperzine, fr_priamide, fr_prisulfonamd,
 fr_pyridine, fr_quatN, fr_sulfide, fr_sulfonamd, fr_sulfone, fr_term_acetylene, fr_tetrazole, fr_thiazole, fr_thiocyan,
 fr_thiophene, fr_unbrch_alkane, fr_urea)


def calc_rdkit_descriptors(mol):
    AllChem.ComputeGasteigerCharges(mol)
    if np.isnan(MinPartialCharge(mol)):
        print("NaN!", Chem.MolToSmiles(mol))
        descriptors = None
    else:
        descriptors = [
                BalabanJ(mol) , # 0
                0.0001*BertzCT(mol) , # 1
                0.1*Chi0(mol) , # 2
                0.1*Chi0n(mol) , # 3
                0.1*Chi0v(mol) , # 4
                0.1*Chi1(mol) , # 5
                0.1*Chi1n(mol) , # 6
                0.1*Chi1v(mol) , # 7
                0.1*Chi2n(mol) , # 8
                0.1*Chi2v(mol) , # 9
                0.1*Chi3n(mol) , # 10
                0.1*Chi3v(mol) , # 11
                0.1*Chi4n(mol) , # 12
                0.1*Chi4v(mol) , # 13
                0.01*EState_VSA1(mol) , # 14
                0.01*EState_VSA10(mol) , #15
                0.01*EState_VSA11(mol) , #16
                0.01*EState_VSA2(mol) , #17
                0.01*EState_VSA3(mol) , #18
                0.01*EState_VSA4(mol) ,
                0.01*EState_VSA5(mol) , #20
                0.01*EState_VSA6(mol) ,
                0.01*EState_VSA7(mol) ,
                0.01*EState_VSA8(mol) ,
                0.01*EState_VSA9(mol) ,
                0.001*ExactMolWt(mol) , #25
                FractionCSP3(mol) ,
                HallKierAlpha(mol) ,
                0.01*HeavyAtomCount(mol) ,
                0.001*HeavyAtomMolWt(mol) ,
                0.1*Kappa1(mol) , #30
                0.1*Kappa2(mol) ,
                0.001*Kappa3(mol) ,
                0.01*LabuteASA(mol) ,
                0.1*MaxAbsEStateIndex(mol) ,
                MaxAbsPartialCharge(mol) , #35
                0.1*MaxEStateIndex(mol) ,
                MaxPartialCharge(mol) , #37
                MinAbsEStateIndex(mol) ,
                MinAbsPartialCharge(mol) , #39
                MinEStateIndex(mol) , #40
                MinPartialCharge(mol) ,
                0.1*MolLogP(mol) ,
                0.01*MolMR(mol) ,
                0.001*MolWt(mol) ,
                0.1*NHOHCount(mol) , #45
                0.1*NOCount(mol) ,
                NumAliphaticCarbocycles(mol) ,
                NumAliphaticHeterocycles(mol) ,
                0.1*NumAliphaticRings(mol) ,
                NumAromaticCarbocycles(mol) , #50
                NumAromaticHeterocycles(mol) ,
                NumAromaticRings(mol) ,
                0.1*NumHAcceptors(mol) ,
                0.1*NumHDonors(mol) ,
                0.1*NumHeteroatoms(mol) , #55
                NumRadicalElectrons(mol) ,
                0.1*NumRotatableBonds(mol) ,
                NumSaturatedCarbocycles(mol) ,
                NumSaturatedHeterocycles(mol) ,
                0.1*NumSaturatedRings(mol) , #60
                0.01*NumValenceElectrons(mol) ,
                0.01*PEOE_VSA1(mol) ,
                0.01*PEOE_VSA10(mol) ,
                0.01*PEOE_VSA11(mol) ,
                0.01*PEOE_VSA12(mol) , #65
                0.01*PEOE_VSA13(mol) ,
                0.01*PEOE_VSA14(mol) ,
                0.01*PEOE_VSA2(mol) ,
                0.01*PEOE_VSA3(mol) ,
                0.01*PEOE_VSA4(mol) , #70
                0.01*PEOE_VSA5(mol) ,
                0.01*PEOE_VSA6(mol) ,
                0.01*PEOE_VSA7(mol) ,
                0.01*PEOE_VSA8(mol) ,
                0.01*PEOE_VSA9(mol) , # 75
                0.1*RingCount(mol) ,
                0.01*SMR_VSA1(mol) ,
                0.01*SMR_VSA10(mol) ,
                0.01*SMR_VSA2(mol) ,
                0.01*SMR_VSA3(mol) , # 80
                0.01*SMR_VSA4(mol) ,
                0.01*SMR_VSA5(mol) ,
                0.01*SMR_VSA6(mol) ,
                0.01*SMR_VSA7(mol) ,
                0.01*SMR_VSA8(mol) , #85
                0.01*SMR_VSA9(mol) ,
                0.01*SlogP_VSA1(mol) ,
                0.01*SlogP_VSA10(mol) ,
                0.01*SlogP_VSA11(mol) ,
                0.01*SlogP_VSA12(mol) ,#90
                0.01*SlogP_VSA2(mol) , #91
                0.01*SlogP_VSA3(mol) ,
                0.01*SlogP_VSA4(mol) ,
                0.01*SlogP_VSA5(mol) ,
                0.01*SlogP_VSA6(mol) , #95
                0.01*SlogP_VSA7(mol) ,
                0.01*SlogP_VSA8(mol) , #97
                0.01*SlogP_VSA9(mol) ,
                0.01*TPSA(mol) ,
                0.01*VSA_EState1(mol) ,#100
                0.01*VSA_EState10(mol) ,
                0.01*VSA_EState2(mol) ,
                0.01*VSA_EState3(mol) ,
                0.01*VSA_EState4(mol) ,
                0.01*VSA_EState5(mol) , #105
                0.01*VSA_EState6(mol) ,
                0.01*VSA_EState7(mol) ,
                0.01*VSA_EState8(mol) ,
                0.01*VSA_EState9(mol) ,
                fr_Al_COO(mol) , #110
                fr_Al_OH(mol) ,
                fr_Al_OH_noTert(mol) ,
                fr_ArN(mol) ,
                fr_Ar_COO(mol) ,
                fr_Ar_N(mol) ,
                fr_Ar_NH(mol) ,
                fr_Ar_OH(mol) ,
                fr_COO(mol) ,
                fr_COO2(mol) ,
                fr_C_O(mol) ,
                fr_C_O_noCOO(mol) ,
                fr_C_S(mol) ,
                fr_HOCCN(mol) ,
                fr_Imine(mol) ,
                fr_NH0(mol) ,
                fr_NH1(mol) ,
                fr_NH2(mol) ,
                fr_N_O(mol) ,
                fr_Ndealkylation1(mol) ,
                fr_Ndealkylation2(mol) ,
                fr_Nhpyrrole(mol) ,
                fr_SH(mol) ,
                fr_aldehyde(mol) ,
                fr_alkyl_carbamate(mol) ,
                fr_alkyl_halide(mol) ,
                fr_allylic_oxid(mol) ,
                fr_amide(mol) ,
                fr_amidine(mol) ,
                fr_aniline(mol) ,
                fr_aryl_methyl(mol) ,
                fr_azide(mol) ,
                fr_azo(mol) ,
                fr_barbitur(mol) ,
                fr_benzene(mol) ,
                fr_benzodiazepine(mol) ,
                fr_bicyclic(mol) ,
                fr_diazo(mol) ,
                fr_dihydropyridine(mol) ,
                fr_epoxide(mol) ,
                fr_ester(mol) ,
                fr_ether(mol) ,
                fr_furan(mol) ,
                fr_guanido(mol) ,
                fr_halogen(mol) ,
                fr_hdrzine(mol) ,
                fr_hdrzone(mol) ,
                fr_imidazole(mol) ,
                fr_imide(mol) ,
                fr_isocyan(mol) ,
                fr_isothiocyan(mol) ,
                fr_ketone(mol) ,
                fr_ketone_Topliss(mol) ,
                fr_lactam(mol) ,
                fr_lactone(mol) ,
                fr_methoxy(mol) ,
                fr_morpholine(mol) ,
                fr_nitrile(mol) ,
                fr_nitro(mol) ,
                fr_nitro_arom(mol) ,
                fr_nitro_arom_nonortho(mol) ,
                fr_nitroso(mol) ,
                fr_oxazole(mol) ,
                fr_oxime(mol) ,
                fr_para_hydroxylation(mol) ,
                fr_phenol(mol) ,
                fr_phenol_noOrthoHbond(mol) ,
                fr_phos_acid(mol) ,
                fr_phos_ester(mol) ,
                fr_piperdine(mol) ,
                fr_piperzine(mol) ,
                fr_priamide(mol) ,
                fr_prisulfonamd(mol) ,
                fr_pyridine(mol) ,
                fr_quatN(mol) ,
                fr_sulfide(mol) ,
                fr_sulfonamd(mol) ,
                fr_sulfone(mol) ,
                fr_term_acetylene(mol) ,
                fr_tetrazole(mol) ,
                fr_thiazole(mol) ,
                fr_thiocyan(mol) ,
                fr_thiophene(mol),
                fr_unbrch_alkane(mol) ,
                fr_urea(mol) , #rdkit properties # 196
                ]
    return descriptors


# In[ ]:


desc_list = []
for m in tqdm(mols):
  desc = calc_rdkit_descriptors(m)
  desc_list.append(desc)


# In[ ]:


desc_df = pd.DataFrame(desc_list)
desc_df


# In[ ]:


from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from rdkit.Chem.Descriptors import FpDensityMorgan1,FpDensityMorgan2,FpDensityMorgan3    ,ExactMolWt,HeavyAtomMolWt,MaxAbsPartialCharge,MaxPartialCharge,MinAbsPartialCharge    ,MinPartialCharge,NumRadicalElectrons,NumValenceElectrons
from rdkit.Chem.rdMolDescriptors import CalcFractionCSP3, CalcKappa1, CalcKappa2, CalcKappa3    ,CalcLabuteASA,CalcNumAliphaticCarbocycles,CalcNumAliphaticHeterocycles    ,CalcNumAliphaticRings,CalcNumAmideBonds,CalcNumAromaticCarbocycles    ,CalcNumAromaticHeterocycles,CalcNumAromaticRings,CalcNumAtomStereoCenters    ,CalcNumBridgeheadAtoms,CalcNumHBA,CalcNumHBD,CalcNumHeteroatoms,CalcNumHeterocycles    ,CalcNumLipinskiHBA,CalcNumLipinskiHBD,CalcNumRings,CalcNumRotatableBonds    ,CalcNumSaturatedCarbocycles,CalcNumSaturatedHeterocycles,CalcNumSaturatedRings    ,CalcNumSpiroAtoms,CalcNumUnspecifiedAtomStereoCenters,CalcTPSA
from subprocess import Popen, PIPE


# In[ ]:


def calc_rdkit_descriptors15(mol):
    AllChem.ComputeGasteigerCharges(mol)
    if np.isnan(MinPartialCharge(mol)):
        print("NaN!", Chem.MolToSmiles(mol))
        descriptors = None
    else:
        descriptors = [
                NumRadicalElectrons(mol) , # 0
                CalcFractionCSP3(mol) , # 1
                CalcNumAliphaticCarbocycles(mol) , # 2
                CalcNumAliphaticHeterocycles(mol) , # 3
                CalcNumAmideBonds(mol) , # 4
                CalcNumAromaticCarbocycles(mol) , # 5
                CalcNumAromaticHeterocycles(mol) , # 6
                CalcNumHBA(mol) , # 7
                CalcNumHBD(mol) , # 8
                CalcNumHeteroatoms(mol) , # 9
                CalcNumHeterocycles(mol) , # 10
                CalcNumLipinskiHBA(mol) , # 11
                CalcNumLipinskiHBD(mol) , # 12
                CalcNumRings(mol) , # 13
                CalcNumRotatableBonds(mol) , # 14
                ]
    return descriptors


# In[ ]:


desc15_list = []
for m in tqdm(mols):
  desc15 = calc_rdkit_descriptors15(m)
  desc15_list.append(desc15)


# In[ ]:


desc15_df = pd.DataFrame(desc15_list)
desc15_df


# In[ ]:


from rdkit import Chem
from mordred import Calculator, descriptors
mols = [Chem.MolFromSmiles(x) for x in df["SMILES"]]
# create descriptor calculator with all descriptors
calc = Calculator(descriptors, ignore_3D=True)

# Commented out IPython magic to ensure Python compatibility.
mordred_df = calc.pandas(mols)

mordred_df




# In[ ]:


object_dtype = [ idx for idx, i in enumerate(mordred_df.dtypes) if i == 'object' ]
new_mordred_df = mordred_df.drop(mordred_df.columns[object_dtype], axis=1)

new_mordred_df


# In[17]:


##model 학습##

X = fp_df #maccs_df #desc_df  #fp_df #new_mordred_df
y = score_df #AD_GPU Top score / best_score



# In[18]:


import sklearn.model_selection
import sklearn.ensemble
import optuna
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)


# In[ ]:


def objective(trial):
    n_estimators = trial.suggest_int('n_estimators', 2, 200)
    max_depth = int(trial.suggest_float('max_depth', 1, 32, log=True))
    min_samples_split = trial.suggest_int('min_samples_split', 2, 15)
    min_samples_leaf = trial.suggest_int('min_samples_leaf', 1, 10)
    max_features = trial.suggest_categorical('max_features', ['auto', 'sqrt', 'log2'])
    clf = RandomForestRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        min_samples_leaf=min_samples_leaf,
        max_features=max_features
    )
    return cross_val_score(clf, X_train, y_train,
                           n_jobs=-1, cv=3,
                           scoring='neg_mean_squared_error').mean()
study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=100)
trial = study.best_trial
print('MSE: {}'.format(trial.value))
print("Best hyperparameters: {}".format(trial.params))


# In[ ]:


from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
# Optuna를 통해 얻은 최적의 하이퍼파라미터
best_params = trial.params
# 모델 생성
model = RandomForestRegressor(
    n_estimators=best_params['n_estimators'],
    max_depth=best_params['max_depth'],
    min_samples_split=best_params['min_samples_split'],
    min_samples_leaf=best_params['min_samples_leaf'],
    max_features=best_params['max_features'],
    n_jobs=2,
    random_state=42)
# 모델 학습
model.fit(X_train, y_train)
# 예측 및 성능 평가
y_pred = model.predict(X_test)


# In[ ]:


##no optimization

my_model = RandomForestRegressor(n_jobs=2, random_state=42)
my_model.fit(X_train, y_train)

y_pred = my_model.predict(X_test)


# In[ ]:


plt.scatter(y_test, y_pred, s=3, alpha=0.8)
plt.xlabel("True docking score", fontsize='xx-large')
plt.ylabel("Predicted docking score", fontsize='xx-large')
plt.grid()
plt.plot(range(-15, 0), range(-15, 0), "r--", linewidth=0.8)
plt.xlim(-15, 0)
plt.ylim(-15, 0)
plt.axis('square')
plt.savefig("50kplot_ecfp2_2048_optim.pdf")


# In[ ]:



from sklearn.metrics import mean_squared_error
mse1 = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error of this model is {mse1:.3f} (docking score unit)")

np.corrcoef(y_test, y_pred)

print(f"Pearson's correlation coefficient of our model is {np.corrcoef(y_test,y_pred)[0,1]:.3f}")


# In[ ]:
