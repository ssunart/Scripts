from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
import matplotlib.pyplot as plt
from PIL import Image
import io

# SMILES 문자열로부터 분자 생성
smiles1 = 'CC(=O)C1=CC=C(OC2=CC(=C(O)C=C2)C(C)(C)C)C=C1'  # 예시 SMILES 문자열
smiles2 = 'Cc3nc2ccc(Oc1ccc(O)c(C(C)(C)C)c1)cc2n3C'       # 다른 예시 SMILES 문자열

mol1 = Chem.MolFromSmiles(smiles1)
mol2 = Chem.MolFromSmiles(smiles2)

# 분자의 지문 생성
fp1 = AllChem.GetMorganFingerprint(mol1, 2)
fp2 = AllChem.GetMorganFingerprint(mol2, 2)

# Tanimoto 유사도 계산
similarity = TanimotoSimilarity(fp1, fp2)
print(f'Tanimoto Similarity: {similarity}')

# 공통 구조 탐색 함수 정의
def find_common_substructures(mol1, mol2):
    mol1_fragments = []
    mol2_fragments = []

    for atom1 in range(mol1.GetNumAtoms()):
        env1 = Chem.FindAtomEnvironmentOfRadiusN(mol1, 2, atom1)
        if env1:
            submol1 = Chem.PathToSubmol(mol1, env1)
            mol1_fragments.append(submol1)

    for atom2 in range(mol2.GetNumAtoms()):
        env2 = Chem.FindAtomEnvironmentOfRadiusN(mol2, 2, atom2)
        if env2:
            submol2 = Chem.PathToSubmol(mol2, env2)
            mol2_fragments.append(submol2)

    common_substructures = []
    for frag1 in mol1_fragments:
        for frag2 in mol2_fragments:
            if frag1.HasSubstructMatch(frag2):
                common_substructures.append(frag1)
    return common_substructures

# 공통 구조 찾기
common_substructures = find_common_substructures(mol1, mol2)

# 공통 구조 정보 출력
print(f"Number of Common Substructures: {len(common_substructures)}")

# 2D 구조 시각화 및 하이라이트 이미지 생성
def get_highlighted_image(mol, substructures):
    highlight_atoms = set()
    for submol in substructures:
        for atom in submol.GetAtoms():
            highlight_atoms.add(atom.GetIdx())
    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
    drawer.DrawMolecule(mol, highlightAtoms=list(highlight_atoms))
    drawer.FinishDrawing()
    image = Image.open(io.BytesIO(drawer.GetDrawingText()))
    return image

# 하이라이트된 이미지 생성
mol1_image_highlight = get_highlighted_image(mol1, common_substructures)
mol2_image_highlight = get_highlighted_image(mol2, common_substructures)

# Matplotlib으로 시각화
fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# 원본 분자 구조
axes[0, 0].imshow(Draw.MolToImage(mol1, size=(300, 300)))
axes[0, 0].axis('off')
axes[0, 0].set_title('Molecule 1')

axes[0, 1].imshow(Draw.MolToImage(mol2, size=(300, 300)))
axes[0, 1].axis('off')
axes[0, 1].set_title('Molecule 2')

# SMARTS 하이라이트 구조
axes[1, 0].imshow(mol1_image_highlight)
axes[1, 0].axis('off')
axes[1, 0].set_title('Molecule 1 Highlighted')

axes[1, 1].imshow(mol2_image_highlight)
axes[1, 1].axis('off')
axes[1, 1].set_title('Molecule 2 Highlighted')

plt.suptitle(f'Tanimoto Similarity: {similarity:.2f}', fontsize=14)
plt.tight_layout()
plt.show()
