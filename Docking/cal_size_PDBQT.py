###PDBQT파일들을 읽고 리간드 크기 분석이나 리간드 좌표의 제한 조건 확인###

import os
import sys
from glob import glob

dic = '/home/soowon/PL_docking/Proj01_DockingBenchmark/Benchmark_input_set/'
files1 = [ file for file in glob(f'{dic}/*/*/*') if '.pdbqt' in file ]
files2 = [ file for file in files1 if 'receptor' in file ]

files = [file for file in files1 if file not in files2]

X, Y, Z = 0,0,0
for file in files:
    cord = [[],[],[]]
    for line in open(file).readlines():
        if line.startswith('ATOM'):
            x, y, z = line[30:38], line[38:46], line[46:54]
            cord[0].append(float(x))
            cord[1].append(float(y))
            cord[2].append(float(z))
    val_x = max(cord[0]) - min(cord[0])
    val_y = max(cord[1]) - min(cord[1])   
    val_z = max(cord[2]) - min(cord[2])
 
    if val_x > 18:
        print(file)
    if X < val_x:
        X = val_x
    if Y < val_y:
        Y = val_y
    if Z < val_z:
        Z = val_z

print(X,Y,Z)

##출력 예시
#특정 파일의 x축 범위가 18 이상이면 파일 경로가 출력됩니다.
#모든 파일 중 x, y, z의 최대 범위가 각각 출력됩니다
