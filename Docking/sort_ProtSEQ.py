###sort Protein sequence from PDB file
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
import csv
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)



# 체인을 선택하는 클래스 정의
class SelectChain(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        if chain.id == self.chain_id:
            return 1
        else:
            return 0

    def accept_residue(self, residue):
        return is_aa(residue)

# PDB 파일 파싱
parser = PDBParser()
structure = parser.get_structure("1ONP", "1ONP.pdb")

# 체인별로 PDB 파일 저장
io = PDBIO()
io.set_structure(structure)

# CSV 파일에 저장할 데이터 초기화
csv_data = []

for model in structure:
    for chain in model:
        io.save(f"chain_{chain.id}.pdb", SelectChain(chain.id))

        # 아미노산 시퀀스 추출
        all_residues = [res for res in chain.get_residues() if is_aa(res)]
        sequence = ""
        for res in all_residues:
            three_letter_code = res.get_resname()
            one_letter_code = protein_letters_3to1.get(three_letter_code, 'X')  # 3-레터 코드를 1-레터 코드로 변환      
            sequence += one_letter_code
        
        # CSV 데이터에 추가
        csv_data.append({"pdbid": "RECEPTOR", "prot_seq": sequence})

# CSV 파일에 저장
with open('protein_sequences.csv', 'w', newline='') as csvfile:
    fieldnames = ['pdbid', 'prot_seq']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for row in csv_data:
        writer.writerow(row)
