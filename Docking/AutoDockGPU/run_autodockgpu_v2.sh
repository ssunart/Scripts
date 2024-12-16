#!/bin/bash
# #SBATCH --partition=g3090_veryshort # gpu[2-4]
#SBATCH --partition=a5000_short # gpu 5
# #SBATCH --partition=quadro_veryshort  # gpu 1
#SBATCH --ntasks=4
#SBATCH --mem=0
#SBATCH --gres=gpu:1
#SBATCH --job-name=Z_3EYI_mn_p19_2
#SBATCH --output=Z_3EYI_mn_p19_2.out


###PART 수정 필요###
##DB PATH chk
module load compiler/gcc-9.2.0 cuda/11.3
TARGET_PATH=/home/soowon/project/06_ZBP1/docking/3EYI/monomer
Lig_PATH=/home/soowon/Database/Enamine/SCR_collection/mk_prep/part19

cp /home/soowon/project/06_ZBP1/docking/scripts/sort_gpuscore.py /home/soowon/project/06_ZBP1/docking/3EYI/monomer/part19_gpuR/dlgXml
    
#시작시간 기록
start_time=$(date +%s)
echo "Start processing at $(date)"

###AutoDock-GPU 수행###
##PART chk##
while read line; do
    echo "Checking file for ${line}"
    # '"$line"'.dlg' 파일이 없으면 다음 명령들을 진행해라.
    if [[ ! -f "part19_gpuR/dlgXml/${line}.dlg" ]]; then
        echo "File ${line}.dlg does not exist, proceeding with autodock."
        autodock_gpu_128wi --lfile "${Lig_PATH}"/"$line".pdbqt \
                           --ffile "${TARGET_PATH}"/3EYI_monomerPrep.maps.fld  \
                           --nrun 50 \
                           --resnam "$line" --xmloutput 0
        
        echo "$line".dlg finished
        grep '^DOCKED' "$line".dlg | cut -c9- > part19_gpuR/pdbqt/"$line"_gpuR.pdbqt
        mv "$line".dlg part19_gpuR/dlgXml/

    fi
done < lig_list_p19.txt


####
cd part19_gpuR/dlgXml
python sort_gpuscore.py
mv top_scores.csv /home/soowon/project/06_ZBP1/docking/3EYI/monomer/part19_gpuR/csv/
#####


# 작업 종료 시간 기록

end_time=$(date +%s)
echo finished at $(date)
                    
# 소요 시간 계산
elapsed_time=$((end_time - start_time))
echo "Time taken for part19: $elapsed_time seconds"

echo "Processing of part19 completed in $elapsed_time seconds" >> processing_times_part19.log
