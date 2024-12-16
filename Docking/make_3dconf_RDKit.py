import os
import csv
from rdkit import Chem
from rdkit.Chem import AllChem

# 1. 현재 디렉토리에 있는 모든 sdf file을 읽어 smiles 정보를 추출한다.
sdf_files = [f for f in os.listdir() if f.endswith('.sdf')]

all_smiles = []

for sdf_file in sdf_files:
    # SDF 파일 읽기
    suppl = Chem.SDMolSupplier(sdf_file)
    
    for mol in suppl:
        # sdf 파일에서 smiles 정보 추출
        smiles = Chem.MolToSmiles(mol)
        all_smiles.append((sdf_file[:-4], smiles))
        
        # 2-2. RDkit을 이용해 새로운 3d ligand sdf file을 생성한다.
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        writer_3d = Chem.SDWriter(f"{sdf_file[:-4]}_3D.sdf")
        writer_3d.write(mol_3d)
        writer_3d.close()
        # all_smiles를 리간드 이름 순으로 정렬
all_smiles = sorted(all_smiles, key=lambda x: int(x[0].split("_")[-1]))
# 2-1. 읽어들인 smiles정보를 해당 리간드 파일 이름과 함께 csv file에 저장한다.
with open('ligands_smiles.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Ligand Name", "SMILES"])
    csvwriter.writerows(all_smiles)
