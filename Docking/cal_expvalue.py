import os, sys, math, glob
import csv

# CSV 파일에 저장할 헤더
header = ['ligand_ID', 'Type', 'Value(nM)', 'calculated']

# 결과를 저장할 리스트
results = []

# 현재 디렉토리의 모든 .sdf 파일을 읽음
for sdf_file in glob.glob('*.sdf'):
    # 파일명에서 확장자를 제외한 부분을 ligand_id로 사용
    ligand_id = os.path.splitext(sdf_file)[0]
    
    # .sdf 파일 읽기
    sortIc50 = open(sdf_file, 'r').readlines()
    
    types = ['Ki (nM)', 'Kd (nM)', 'IC50 (nM)', 'EC50 (nM)']
    
    for t in types:
        a = []
        b = []
        
        for idx, line in enumerate(sortIc50):
            if line.startswith('> <BindingDB MonomerID>'):
                start1 = idx + 1
                a.append(start1)
            if line.startswith(f'> <{t}>'):
                start2 = idx + 1
                b.append(start2)
        
        for idx, i in enumerate(a):
            score = sortIc50[b[idx]].strip()
             # 부등호 제거
            score = score.replace('<', '').replace('>', '')
            
            if score != '':
                start3 = float(score)
                r_score = math.log(start3 * 10**(-9))*(300)*(1.987)*(0.001)
                
                # 결과 저장
                results.append([ligand_id, t, f'{score}(nM)', r_score])

# Ligand_ID로 결과 정렬
results = sorted(results, key=lambda x: (int(''.join(filter(str.isdigit, x[0])) if any(char.isdigit() for char in x[0]) else 0), x[0]))

# 결과를 하나의 CSV 파일에 저장
with open('combined.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # 헤더 쓰기
    csvwriter.writerow(header)
    
    # 결과 쓰기
    csvwriter.writerows(results)
