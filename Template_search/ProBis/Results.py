###ProBiS 탐색 결과 추출
#JSON 데이터를 읽고 특정 기준(예: alignment_score > 10)을 만족하는 데이터를 추출하여 처리한 후, 결과를 CSV 파일로 저장

import json
import csv
with open('example.json', 'r') as f:
    data = json.load(f)

# 높은 alignment 점수를 가진 단백질을 찾는 코드

high_score_proteins = []

# 각 단백질과 그에 대한 alignment 점수를 순회
for protein in data:
    pdb_id = protein['pdb_id']
    chain_id = protein['chain_id']
    alignments = protein['alignment']
    
    for alignment in alignments:
        score = alignment['scores']['alignment_score']
        
        # 예를 들어, alignment 점수가 10 이상인 경우에만 저장
        if score > 10:
            high_score_proteins.append({
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'alignment_score': score
            })

# 결과 출력 (또는 다른 분석을 위해 사용)
print(high_score_proteins)

# CSV 파일로 저장
csv_columns = ['pdb_id', 'chain_id', 'alignment_score']
csv_file = "high_score_proteins.csv"
try:
    with open(csv_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for data in high_score_proteins:
            writer.writerow(data)
except IOError:
    print("I/O error")
