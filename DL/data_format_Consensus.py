###
import csv

# 데이터 수집 예시 (실제 데이터로 대체하세요)
data = [
    {'Set No.': 'set001', 'Protein seq': 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAN', 'Prot family': 'Kinase', 'Lig No.': 'lig001', 'Lig smiles': 'CCOCCOC(=O)CC', 'Exp_value': 7.5},
    {'Set No.': 'set002', 'Protein seq': 'MKGLDIQKRLVEEFKREKA', 'Prot family': 'Phosphatase', 'Lig No.': 'lig002', 'Lig smiles': 'CCCCCC(=O)OC1=CC=CC=C1', 'Exp_value': 5.2},
    # ... 추가 데이터
]

# 결과를 CSV 파일에 저장
with open('protein_ligand_data.csv', 'w', newline='') as csvfile:
    fieldnames = ['Set No.', 'MainPDB','Protein seq', 'Prot family', 'Lig No.', 'Lig smiles', 'Exp_value']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()  # 헤더 작성
    for row in data:
        writer.writerow(row)  # 각 행의 데이터 작성
