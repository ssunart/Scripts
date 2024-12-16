#!/bin/bash



prepare_receptor -r ${TARGET}.pdb -A 'checkhydrogens' -U 'nphs_waters' -e -o ${TARGET}_prep.pdbqt

mk_prepare_ligand.py -i ${TARGET}_Nlig.mol2 -o ${TARGET}_Nlig.pdbqt


### Column 정리 주의 ###
pythonsh prepare_gpf4.py -l ${TARGET}_Nlig.pdbqt -r ${TARGET}_prep.pdbqt \
                         -p ligand_types='HD,Br,Cl,A,C,N,P,NA,OA,F,S,SA,I' \
                         -p npts='50, 50, 50' -p spacing='0.4' -y

autogrid4 -p ${TARGET}_prep.gpf -l ${TARGET}_prep.glg

while read line
do
    autodock_gpu_128wi --lfile ligands/"$line".pdbqt  --ffile ${TARGET}_prep.maps.fld \
                       --nrun 100 --resnam "$line"_gpuR
    grep '^DOCKED' "$line"_gpuR.dlg | cut -c9- > "$line"_gpuR.pdbqt
    mv *"$line"_gpuR* gpuR/
done < ligands/sdf_list
