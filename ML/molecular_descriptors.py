##RDKit 을 이용해서 분자 읽어들이기.
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem.Draw import IPythonConsole # Jupyter notebook에 분자가 바로 표현되게 하려면 이 줄이 필요하다. 

### 분자를 SMILES로 부터 읽어들이기.
m=Chem.MolFromSmiles('C[C@H](O)c1ccccc1')

###SMILES에서 chirality 표시하기.
m2=Chem.MolFromSmiles('CC(O)c1ccccc1')
Chem.MolToSmiles(m2)

#Aromatic 원자의 경우, 소문자로 표현하는 것이 기본이다.
#그러나 필요에 따라서 Kekule form으로 표현하도록 할 수 있다.
Chem.MolToSmiles(m2, kekuleSmiles=True)
m3 = Chem.MolFromSmiles('C1=CC=CN=C1')
Chem.MolToSmiles(m3)
Chem.Kekulize(m3) # aromatic bond가 single-double-single-double-... 이런식으로 표현되도록 한다. 
Chem.MolToSmiles(m3, kekuleSmiles=True)

#SMILES로 출력하기.
Chem.MolToSmiles(m)
Chem.MolToSmiles(m2)

#Kekulization
##single-double-single-... 이런 형태로 SMILES가 출력하도록 함.
Chem.MolToSmiles(m, kekuleSmiles=True, canonical=True) # Canonical SMILES를 출력한다.
Chem.MolToSmiles(m3, kekuleSmiles=True)
Chem.MolToSmiles(m3, kekuleSmiles=False)

#isomericSmilles = True 이면 chirality 정보를 출력한다.
Chem.MolToSmiles(m, kekuleSmiles=True, isomericSmiles=True) 
#isomericSmilles = False 이면 chiraity 정보 무시.
Chem.MolToSmiles(m, kekuleSmiles=True, isomericSmiles=False) 

##분자를 Mol format (SDF format)으로 출력하기.
print(Chem.MolToMolBlock(m)) # rdkit mol instance를 Mol-format으로 출력. 
print(Chem.MolToMolBlock(m2)) # rdkit mol instance를 Mol-format으로 출력. 
Chem.MolToMolFile(m3, "m3.mol")

##PDB 파일로 저장하기
print(Chem.MolToPDBBlock(m3))

from rdkit.Chem import AllChem # AllChem module 읽어들이기.
m3 = Chem.AddHs(m3) # 수소 붙이기
AllChem.EmbedMolecule(m3) # 분자의 3차원 구조 생성
print(Chem.MolToPDBBlock(m3))

Chem.MolToPDBFile(m3, "m3.pdb")

###수소 붙이기
m3_with_H = Chem.AddHs(m3) 
print(Chem.MolToMolBlock(m3_with_H))

##분자가 복잡할 경우, 수소의 좌표가 이상하게 붙는 경우들이 존재한다.
##이럴 경우, 수소의 좌표를 다시 계산해주어야 한다.
##AllChem.Compute2DCoords 함수를 사용하여야 한다.

from rdkit.Chem import AllChem # AllChem module 읽어들이기.
AllChem.Compute2DCoords(m4_with_H) # m4_with_H 그 자체에서 좌표가 계산되어서 저장되었음. 
m4_with_H
print(Chem.MolToMolBlock(m4_with_H))

#3D 좌표에 관해
"""주의: 모든 z 축이 0으로 설정되어 있다.
기본적으로 그림으로 표현된 분자식 기반의 좌표를 저장하고 있기 때문에 2차원으로 표현되어 있음.
3차원 구조를 얻어내기 위해서는 추가적인 계산이 필요함.
하지만 분자가 커질수록 최소 에너지를 가지는 분자의 3차원 구조를 찾는 것이 쉽지 않다.
근사적으로 3^(Number of rotatable bonds) 개수만큼 3D 구조가 가능하다.
그 중에서 에너지가 가장 낮은 것을 찾는 것은 쉽지 않다.
그래서 rdkit에서는 약간의 가정과 근사를 써서 3D 구조를 만들어낸다."""
AllChem.EmbedMolecule(m4_with_H) # 분자 클래스 자체를 변형시킨다. 

#Force-field를 이용한 분자의 3차원 구조 최적화
#분자의 3차원 구조를 한 번 더 최적화 하기 위해서는 아래의 AllChem.MMFFOptimizeMolecule 함수를 사용한다.
AllChem.MMFFOptimizeMolecule(m4_with_H)

#3D 구조 표현에서 2D 구조로 바꾸기
AllChem.Compute2DCoords(m4_with_H)

##Fingerprint 계산
m1 = Chem.MolFromSmiles('Cc1ccccc1')

#Extended Fingerprint (ECFP) 계산
fp1 = AllChem.GetMorganFingerprint(m1, 2) # radius = 2, diameter = 4 이므로 ECFP4 에 해당. 

#Morgan Fingerprint as a bit vector
fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=1024) # radius = 2, bit = 1024. 

#Total polar surface area 계산 (Angstrom^2 단위)
from rdkit.Chem import Descriptors
m = Chem.MolFromSmiles('c1ccccc1C(=O)O')
Descriptors.TPSA(m)

#logP 계산
Descriptors.MolLogP(m)

#부분 전하
m = Chem.MolFromSmiles('c1ccccc1C(=O)O')
AllChem.ComputeGasteigerCharges(m) # 가장 많이 사용되는 부분 전하 계산 방법, GasteigerCharge
m.GetAtomWithIdx(0).GetDoubleProp('_GasteigerCharge')

##각 descriptor의 contribution을 계산 및 시각화
from rdkit.Chem.Draw import SimilarityMaps
mol = Chem.MolFromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')
AllChem.ComputeGasteigerCharges(mol) # 부분전하 계산
contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
##높은 값을 가지는 원자는 빨간색을 띄고, 낮은 값을 가지는 원자는 파란색으로 표현된다.

#logP의 기여도를 시각화 해보자.
from rdkit.Chem import rdMolDescriptors
contribs = rdMolDescriptors._CalcCrippenContribs(mol)
fig = SimilarityMaps.GetSimilarityMapFromWeights(mol,[x for x,y in contribs], colorMap='jet', contourLines=10)


#RDKit을 이용한 195가지의 descriptor 계산

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
        finger = None
    else:
        finger = [
                BalabanJ(mol) , # 0
                0.0001*BertzCT(mol) , # 1
                0.1*Chi0(mol) , # 2
                0.1*Chi0n(mol) , # 3
                0.1*Chi0v(mol) , # 4
                0.1*Chi1(mol) , # 5
                0.1*Chi1n(mol) , # 6
                0.1*Chi1v(mol) , # 7
                0.1*Chi2n(mol) ,
                0.1*Chi2v(mol) ,
                0.1*Chi3n(mol) ,
                0.1*Chi3v(mol) ,
                0.1*Chi4n(mol) ,
                0.1*Chi4v(mol) ,
                0.01*EState_VSA1(mol) ,
                0.01*EState_VSA10(mol) ,
                0.01*EState_VSA11(mol) ,
                0.01*EState_VSA2(mol) ,
                0.01*EState_VSA3(mol) ,
                0.01*EState_VSA4(mol) ,
                0.01*EState_VSA5(mol) ,
                0.01*EState_VSA6(mol) ,
                0.01*EState_VSA7(mol) ,
                0.01*EState_VSA8(mol) ,
                0.01*EState_VSA9(mol) ,
                0.001*ExactMolWt(mol) ,
                FractionCSP3(mol) ,
                HallKierAlpha(mol) ,
                0.01*HeavyAtomCount(mol) ,
                0.001*HeavyAtomMolWt(mol) ,
                # Ipc(mol) ,
                0.1*Kappa1(mol) ,
                0.1*Kappa2(mol) ,
                0.001*Kappa3(mol) ,
                0.01*LabuteASA(mol) ,
                0.1*MaxAbsEStateIndex(mol) ,
                MaxAbsPartialCharge(mol) ,
                0.1*MaxEStateIndex(mol) ,
                MaxPartialCharge(mol) ,
                MinAbsEStateIndex(mol) ,
                MinAbsPartialCharge(mol) ,
                MinEStateIndex(mol) ,
                MinPartialCharge(mol) ,
                0.1*MolLogP(mol) ,
                0.01*MolMR(mol) ,
                0.001*MolWt(mol) ,
                0.1*NHOHCount(mol) ,
                0.1*NOCount(mol) ,
                NumAliphaticCarbocycles(mol) ,
                NumAliphaticHeterocycles(mol) ,
                0.1*NumAliphaticRings(mol) ,
                NumAromaticCarbocycles(mol) ,
                NumAromaticHeterocycles(mol) ,
                NumAromaticRings(mol) ,
                0.1*NumHAcceptors(mol) ,
                0.1*NumHDonors(mol) ,
                0.1*NumHeteroatoms(mol) ,
                NumRadicalElectrons(mol) ,
                0.1*NumRotatableBonds(mol) ,
                NumSaturatedCarbocycles(mol) ,
                NumSaturatedHeterocycles(mol) ,
                0.1*NumSaturatedRings(mol) ,
                0.01*NumValenceElectrons(mol) ,
                0.01*PEOE_VSA1(mol) ,
                0.01*PEOE_VSA10(mol) ,
                0.01*PEOE_VSA11(mol) ,
                0.01*PEOE_VSA12(mol) ,
                0.01*PEOE_VSA13(mol) ,
                0.01*PEOE_VSA14(mol) ,
                0.01*PEOE_VSA2(mol) ,
                0.01*PEOE_VSA3(mol) ,
                0.01*PEOE_VSA4(mol) ,
                0.01*PEOE_VSA5(mol) ,
                0.01*PEOE_VSA6(mol) ,
                0.01*PEOE_VSA7(mol) ,
                0.01*PEOE_VSA8(mol) ,
                0.01*PEOE_VSA9(mol) ,
                0.1*RingCount(mol) ,
                0.01*SMR_VSA1(mol) ,
                0.01*SMR_VSA10(mol) ,
                0.01*SMR_VSA2(mol) ,
                0.01*SMR_VSA3(mol) ,
                0.01*SMR_VSA4(mol) ,
                0.01*SMR_VSA5(mol) ,
                0.01*SMR_VSA6(mol) ,
                0.01*SMR_VSA7(mol) ,
                0.01*SMR_VSA8(mol) ,
                0.01*SMR_VSA9(mol) ,
                0.01*SlogP_VSA1(mol) ,
                0.01*SlogP_VSA10(mol) ,
                0.01*SlogP_VSA11(mol) ,
                0.01*SlogP_VSA12(mol) ,
                0.01*SlogP_VSA2(mol) ,
                0.01*SlogP_VSA3(mol) ,
                0.01*SlogP_VSA4(mol) ,
                0.01*SlogP_VSA5(mol) ,
                0.01*SlogP_VSA6(mol) ,
                0.01*SlogP_VSA7(mol) ,
                0.01*SlogP_VSA8(mol) ,
                0.01*SlogP_VSA9(mol) ,
                0.01*TPSA(mol) ,
                0.01*VSA_EState1(mol) ,
                0.01*VSA_EState10(mol) ,
                0.01*VSA_EState2(mol) ,
                0.01*VSA_EState3(mol) ,
                0.01*VSA_EState4(mol) ,
                0.01*VSA_EState5(mol) ,
                0.01*VSA_EState6(mol) ,
                0.01*VSA_EState7(mol) ,
                0.01*VSA_EState8(mol) ,
                0.01*VSA_EState9(mol) ,
                fr_Al_COO(mol) ,
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
                
    return finger

