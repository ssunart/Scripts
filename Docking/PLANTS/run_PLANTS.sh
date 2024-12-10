###SLURM 작업 스케줄러를 통해 PLANTS 프로그램 실행
#설정 파일(plantsconfig)을 기반으로 다수의 리간드를 스크리닝
#작업 시작 및 종료 시간을 기록하고, 결과 로그를 저장(plants.out)

#!/bin/bash
#SBATCH -J T1214_plants
#SBATCH -p skylake_normal
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -o /home/sangmin/casp16/T1214/6V81/plants.out
start_time=$(date +"%Y-%m-%d %H:%M:%S")
end_time=$(date +"%Y-%m-%d %H:%M:%S")

echo "Job started at $start_time"
PLANTS --mode screen /home/sangmin/casp16/T1214/6V81/plantsconfig
echo "Job ended at $end_time"
