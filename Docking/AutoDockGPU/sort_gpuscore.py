import os
import csv
import re

file_path = os.getcwd()
print(f"Score sorting directory: {file_path}")

def read_dlg_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    start_reading = False
    min_binding_energy = float('inf')
    
    for line in lines:
        if "RMSD TABLE" in line:
            start_reading = True
            continue

        if start_reading:
            try:
                # Binding Energy 값을 찾아냅니다.
                binding_energy = float(line.split()[3])
                min_binding_energy = min(min_binding_energy, binding_energy)
            except (ValueError, IndexError):
                continue
                
    return min_binding_energy

# 현재 디렉토리의 모든 파일을 검색
files = [f for f in os.listdir('.') if os.path.isfile(f)]

# 결과를 저장할 리스트
results = []

# 각 .dlg 파일에 대해 작업을 수행
for file in files:
    if file.endswith('.dlg'):
        min_binding_energy = read_dlg_file(file)
        results.append((file, min_binding_energy))

# 결과를 CSV 파일에 저장
with open('top_scores.csv', 'w', newline='') as csvfile:
    fieldnames = ['Lig No.', 'AD_GPU']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for file, score in results:
        writer.writerow({'Lig No.': file, 'AD_GPU': score})

import csv

# CSV 파일에서 데이터를 읽습니다.
with open('top_scores.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)  # 헤더를 건너뜁니다.
    rows = [row for row in reader]

# 리간드 번호(Lig No.)를 기준으로 정렬합니다.
sorted_rows = sorted(rows, key=lambda x: int(x[0].split('_')[2].replace('gpuR.dlg', '')))

# 정렬된 데이터를 다시 CSV 파일에 저장합니다.
with open('sorted_top_scores.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)  # 헤더를 작성합니다.
    writer.writerows(sorted_rows)  # 정렬된 행을 작성합니다.
