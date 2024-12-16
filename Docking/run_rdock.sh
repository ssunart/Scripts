#!/bin/bash

# prep.prm 파일 복사
cp ../../../process/prep.prm .

# cvt_prot.cxc 파일 생성
filename="cvt_prot.cxc"
content="open ../1ONP_prep.pdbqt
save ./1ONP_prep.mol2
exit"
echo "$content" > $filename
echo "$filename 파일이 생성되었습니다."

# ChimeraX 명령 실행
chimerax --nogui cvt_prot.cxc

# Open Babel을 사용해 mol2 파일을 sdf로 변환
obabel -imol2 ../1ONP_Nligand.mol2 -osdf -O 1ONP_Nligand.sdf

# prep.prm 파일 수정 (현재는 실제로 치환이 발생하지 않음)
sed -i 's/TRGT/1ONP/g' prep.prm

# rbcavity 명령 실행
rbcavity --was -d -r prep.prm

# PDBQT 파일들이 위치한 디렉토리 경로
input_dir="../ligands"
# SDF 파일들이 저장될 디렉토리 경로
output_dir="ligands"

# output_dir이 존재하지 않으면 생성
mkdir -p "$output_dir"

# input_dir 내의 모든 pdbqt 파일에 대해 작업 수행
for pdbqt_file in "$input_dir"/*.pdbqt; do
    # 파일명만 추출 (확장자 제외)
    base_name=$(basename "$pdbqt_file" .pdbqt)

    # obabel 명령 실행하여 sdf 파일로 변환
    obabel -ipdbqt "$pdbqt_file" -o sdf -O "$output_dir/$base_name.sdf"

    # 변환 결과 출력
    if [ $? -eq 0 ]; then
        echo "$pdbqt_file -> $output_dir/$base_name.sdf 변환 성공"
    else
        echo "$pdbqt_file 변환 실패"
    fi
done

# SDF 파일들이 위치한 디렉토리 경로
input_dir="ligands"

# input_dir 내의 모든 sdf 파일에 대해 작업 수행
for sdf_file in "$input_dir"/*.sdf; do
    # 파일명만 추출 (확장자 제외)
    base_name=$(basename "$sdf_file" .sdf)

    # rbdock 명령 실행
    rbdock -i "$sdf_file" -o "${base_name}_rD" -r prep.prm -p dock.prm -n 100

    # sdreport 명령 실행
    sdreport -t "${base_name}_rD.sd" > "rDreport_${base_name}.txt"

    # 결과 출력
    if [ $? -eq 0 ]; then
        echo "$sdf_file -> ${base_name}_rD.sd 및 rDreport_${base_name}.txt 생성 성공"
    else
        echo "$sdf_file 작업 실패"
    fi
done
