#!/bin/bash
#SBATCH --partition=g4090_normal # gpu[2-4]
# #SBATCH --partition=a5000_short # gpu 5
# #SBATCH --partition=quadro_veryshort  # gpu 1
#SBATCH --ntasks=4
#SBATCH --mem=0
#SBATCH --gres=gpu:1
#SBATCH --job-name=Z_3EYI_di_p04
#SBATCH --output=Z_3EYI_di_p04.out


###PART 수정 필요###
##DB PATH chk
module load compiler/gcc-9.2.0 cuda/11.3
TARGET_PATH=/home/soowon/project/06_ZBP1/raw_file/3EYI
Lig_PATH=/home/soowon/Database/Enamine/SCR_collection/mk_prep/part4

cp ../scripts/sort_gpuscore.py .
    
#시작시간 기록
start_time=$(date +%s)
echo "Start processing at $(date)"

###AutoDock-GPU 수행###
##PART chk##
while read line; do
    echo "Checking file for ${line}"
    # '"$line"'.dlg' 파일이 없으면 다음 명령들을 진행해라.
    if [[ ! -f "part4_gpuR/dlgXml/${line}.dlg" ]]; then
        echo "File ${line}.dlg does not exist, proceeding with autodock."
        autodock_gpu_128wi --lfile "${Lig_PATH}"/"$line".pdbqt \
                           --ffile "${TARGET_PATH}"/3EYI_dimerPrep.maps.fld  \
                           --nrun 50 \
                           --resnam "$line"
        
        echo "$line".dlg finished
        grep '^DOCKED' "$line".dlg | cut -c9- > "$line"_gpuR.pdbqt
        mv "$line".dlg part4_gpuR/dlgXml/
        mv "$line".xml part4_gpuR/dlgXml/
        mv "$line"_gpuR.pdbqt part4_gpuR/pdbqt/

    fi
done < lig_list_p04.txt

cd part4_gpuR/dlgXml
cp ../sort_gpuscore.py .
python sort_gpuscore.py
mv top_scores.csv ../csv/

# 작업 종료 시간 기록
end_time=$(date +%s)
echo finished at $(date)
                    
# 소요 시간 계산
elapsed_time=$((end_time - start_time))
echo "Time taken for part4: $elapsed_time seconds"

echo "Processing of part4 completed in $elapsed_time seconds" >> processing_times_part4.log
