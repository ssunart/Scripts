###SLURM 작업 스케줄러를 통해 특정 경로에서 DOCK6를 실행하는 스크립트를 생성하고 제출###

import subprocess as sp
import os


def sph(path):
    batch = f"""#!/bin/bash
#SBATCH --partition=skylake
#SBATCH --qos=normal
#SBATCH --ntasks=2
#SBATCH --job-name=sph
#SBATCH --output sph.out

source /home/byun/apps/dock6/dock6.sh
source /home/soowon/.bashrc

cd {path}
sphgen
boxsphere < sphgen_cluster.in
printf "Y\n 5\n rec.sph\n 1\n rec_box.pdb\n" | showboxY
grid -i grid.in
    """
    fileName = f'{path}/tmp.sh'
    f = open(fileName, 'w')
    f.write(batch)
    f.close()
    sp.call(f'sbatch {fileName}', shell=True)
                    
             
paths = ['/home/soowon/PL_docking/UCSF/set574_2qrg', '/home/soowon/PL_docking/UCSF/set568_2w4i']

for path in paths:
    sph(path)
###출력결과: 각 경로(/home/soowon/PL_docking/UCSF/set574_2qrg, /home/soowon/PL_docking/UCSF/set568_2w4i)에서 tmp.sh라는 Bash 스크립트가 생성
###생성된 스크립트는 SLURM에 의해 제출되고 작업 실행
###작업 결과는 각 경로에서 sph.out 파일에 기록
