### SMILES 문자열로부터 분자를 생성하고, 이를 3D 구조로 변환한 후 SDF 파일로 저장
import os


# csv 파일 읽기
data = pd.read_csv('output.smi')

for index, row in data.iterrows():
    catalog_id = row['CatalogID']
    smiles = row['SMILES']

    output_file = f'{catalog_id}.sdf'

    if os.path.isfile(output_file):
        print(f"File {output_file} already exists, Skipping...")
        continue

            # SMILES 문자열로부터 분자 생성
    molecule = Chem.MolFromSmiles(smiles)

        # 3D 구조 생성
    molecule_with_3d_coordinates = Chem.AddHs(molecule)
    embed_status = AllChem.EmbedMolecule(molecule_with_3d_coordinates)

        # EmbedMolecule의 반환값이 0이면 성공, -1이면 실패
    if embed_status != -1:
        AllChem.MMFFOptimizeMolecule(molecule_with_3d_coordinates)

            # sdf 파일로 저장
        writer = Chem.SDWriter(f'{catalog_id}.sdf')
        writer.write(molecule_with_3d_coordinates)
        writer.close()
    else:
        print(f"Embedding failed for CatalogID {catalog_id}")

#Input: output.smi CSV 파일 (각 행에 CatalogID와 SMILES 열 포함).
#Output: 성공적으로 3D 좌표가 생성된 분자는 각각의 SDF 파일(CatalogID.sdf)로 저장
##이미 존재하는 파일은 건너뛰어 중복 저장 방지
##SMILES로부터 3D 구조 생성이 실패할 경우 메시지 출력
