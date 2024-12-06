###특정 경로의 dlg 파일들을 pdbqt 파일로 변환하는 코드###

#!/bin/bash
target=AD_GPU_R
wdir=/home/soowon/\[PL_Dock\]/PROJECT01/BindingDB_set01_process/1onp_ori/${target}/
#wdir=/home/soowon/

cd $wdir
path=`pwd -P`
echo $path
files=$(find $path/ -name "*.dlg")

for file in $files; do
echo ${file%.*}
grep '^DOCKED' ${file} | cut -c9- > ${file%.*}.pdbqt

done

echo  Finish !
