###probis 도구를 사용해 단백질 표면 추출 및 분석
###Python 스크립트를 호출하여 결과(JSON 파일) 확인

#!/bin/bash
shopt -s expand_aliases
alias probis='/home/soowon/apps/probis/probis'

TARGET=TARGET/
probis -extract -f1 ${TARGET}.pdb -c1 A -srffile ${TARGET}.srf

cd protein/

probis -ncpu 2 -surfdb -local -sfile prot_srf.txt -f1 ../${TARGET}.srf -c1 A -nosql example.nosql

probis -results -f1 ../${TARGET}.pdb -c1 A -nosql example.nosql -json example.json
python check_json.py
