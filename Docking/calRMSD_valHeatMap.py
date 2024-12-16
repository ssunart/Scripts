import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 현재 디렉토리의 모든 하위 디렉토리를 탐색합니다.
current_dir = os.getcwd()

for root, dirs, files in os.walk(current_dir):
    # SDF 파일들을 불러올 디렉토리 경로를 설정합니다.
    sdf_files = [os.path.join(root, f) for f in files if f.endswith('.sdf')]
    
    # SDF 파일이 없는 디렉토리는 건너뜁니다.
    if not sdf_files:
        continue
    
    # 파일 이름에서 'ModelX' 부분만 추출합니다.
    sdf_filenames = [os.path.basename(f).split('_')[1] for f in sdf_files]
    
    # 파일 이름과 경로를 함께 저장하여 정렬
    sdf_files_sorted = sorted(zip(sdf_filenames, sdf_files))
    
    # 정렬된 파일 이름과 경로를 분리
    sdf_filenames_sorted, sdf_files_sorted = zip(*sdf_files_sorted)
    
    # 분자 개수를 가져옵니다.
    num_molecules = len(sdf_files_sorted)
    
    # RMSD 값을 저장할 배열을 초기화합니다.
    rmsd_matrix = np.zeros((num_molecules, num_molecules))
    
    # RMSD 값을 계산합니다.
    for i in range(num_molecules):
        for j in range(i + 1, num_molecules):
            mol1 = sdf_files_sorted[i]
            mol2 = sdf_files_sorted[j]
            result = subprocess.run(['obrms', mol1, mol2], capture_output=True, text=True)
            rmsd = float(result.stdout.strip().split()[-1])
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
    
    # 히트맵을 생성합니다.
    plt.figure(figsize=(12, 10))
    sns.heatmap(rmsd_matrix, annot=True, cmap='coolwarm', xticklabels=sdf_filenames_sorted, yticklabels=sdf_filenames_sorted)
    plt.title(f'RMSD Heatmap for {os.path.basename(root)}')
    plt.xlabel('Molecule Index')
    plt.ylabel('Molecule Index')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    
    # 고해상도로 이미지를 저장합니다.
    heatmap_filename = os.path.join(root, 'rmsd_heatmap_high_res.png')
    plt.savefig(heatmap_filename, dpi=300)
    plt.close()
