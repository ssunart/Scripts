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
