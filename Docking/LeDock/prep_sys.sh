#!/bin/bash
###Check Set No. & pdb PDB ID LIG No. ###
#mkdir set499_result
pdb=1p44_Lig_
PDB=1P44

for num in {0..18}
do
     obabel -i sdf ../Lig_split/set499/$pdb$num.sdf -o mol2 -O set499_result/$pdb$num.mol2
done

cd set499_result

lepro_linux_x86 ../../protein/$PDB.pdb
sed -i "13s/.*/100/g" dock.in
