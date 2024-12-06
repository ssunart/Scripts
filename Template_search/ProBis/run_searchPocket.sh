#!/bin/bash
shopt -s expand_aliases
alias probis='/home/soowon/apps/probis/probis'

cd BSDB_231103

# 상위 디렉토리에서 시작
for dir in */ ; do
    # 디렉토리 이름 추출 (마지막 '/' 제거)
    dir_name=${dir%/}

    # 디렉토리로 이동
    cd "$dir"

    # binding pocket 파일들을 찾아서 처리
    for file in ligand_pocket*.pdb; do
        if [ -f "$file" ]; then
            # 파일 이름에서 chain index 추출
            chain=$(echo $file | sed -e 's/ligand_pocket_[0-9]*_\(.*\).pdb/\1/')

            # probis 명령 실행, 파일 이름에 디렉토리 이름 포함
            probis -extract -f1 "$file" -c1 $chain -srffile "${dir_name}_${file%.*}_pocket.srf"
            mv *_pocket.srf ../../SRF_data/pocket
        fi
    done
    # 원래 디렉토리로 돌아가기
    cd ..
done
