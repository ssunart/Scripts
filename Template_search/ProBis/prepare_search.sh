###BS_DB로부터 srf파일 셋업하기

#!/bin/bash
shopt -s expand_aliases
alias probis='/home/soowon/apps/probis/probis'

cd BSDB_231103

# 상위 디렉토리에서 시작
for dir in */ ; do
    dir_name=${dir%/}
    # 디렉토리로 이동
    cd "$dir"

    # ligand_pocket 파일들을 찾아서 처리
    for file in ligand_pocket*.pdb; do
        if [ -f "$file" ]; then
            # 파일 이름에서 chain index 추출
            chain=$(echo $file | sed -e 's/ligand_pocket_[0-9]*_\(.*\).pdb/\1/')

            # probis 명령 실행
            probis -extract -f1 protein.pdb -c1 $chain -srffile "${dir_name}_${file%.*}_prot.srf"
            mv *_prot.srf ../../SRF_data/protein
        fi
    done

    # 원래 디렉토리로 돌아가기
    cd ..
done

