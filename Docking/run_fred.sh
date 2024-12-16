#!/bin/bash
# #SBATCH --partition=rome_short    # cpu1
#SBATCH --partition=skylake_veryshort # cpu[2-7]
#SBATCH --ntasks=4
#SBATCH --mem=0
#SBATCH --job-name=rFr_BindDB_SNUM
#SBATCH --output=rFr_BindDB_SNUM.out



getstructure "VAR_A"
spruce -in "VAR_A".pdb
receptorindu -in *"VAR_B"_A*.oedu -out rec_"VAR_A".oedu

mkdir -p ligands
cp ../ligands/*_Lig*.sdf ligands/
tar -zcvf combined.sdf.gz ligands

fred -mpi_np 4 -receptor rec_"VAR_A".oedu -dbase combined.sdf.gz -num_poses 100 -prefix "VAR_A"_fred
