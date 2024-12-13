#!/bin/bash

shopt -s expand_aliases
alias probis='/home/soowon/apps/probis/probis'

TARGET="/home/soowon/project/05_CASP16/proBis/testDB/EX/8wlb_test"
FULL_NAME=${TARGET##*/}    # '/' 뒤의 문자열 추출: '5sst_test'
DESIRED_NAME=${FULL_NAME%_*}  # '_' 앞 부분만 추출: '5sst'

probis -extract -f1 ${TARGET}.pdb -c1 A -srffile ${TARGET}.srf
cd SRF_data/pocket/
rm *.json
rm *cons.pdb
rm example.nosql


start_time=$(date +%s)
echo "Start processing at $(date)"


probis -ncpu 4 -surfdb -local -sfile ../pock_srf.txt -f1 ${TARGET}.srf -c1 A -nosql example.nosql

probis -results -f1 ${TARGET}.pdb -c1 A -nosql example.nosql -json example.json
python check_json.py

mkdir -p ../testR/${DESIRED_NAME}
mv high_score_proteins.csv ../testR/${DESIRED_NAME}/highscore_pock.csv


end_time=$(date +%s)
echo finished at $(date)
