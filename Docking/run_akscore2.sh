#!/bin/bash
#SBATCH --job-name=AK2_4cff
#SBATCH --ntasks=4
#SBATCH --partition=mix_veryshort
#SBATCH --output=AK2_C1V_7e2e.log
#SBATCH --gres=gpu:1

echo $(date)

module load cuda/11.8

python=/appl/anaconda3/envs/VS_env/bin/python  # python for Akscore2 & RTMScore

Akscore2=/appl/git-repo/AKScore2/run.py     # running Akscore2
RTMScore=/appl/git-repo/RTMScore/rtmscore.py         # running RTMScore
 
mol2pdb=/appl/git-repo/AKScore2/convert_mol2pdb.py   # Making input ligand PDB for Akscore

protein_pdb='./7e2e_prep.pdb'     # INPUT protein pdb for Akscore2 & RTMScore
#PATH_pdb='/home/soowon/project/07_PPAR/rescoring/7e2e_PPAR/dnv_prep'

START_TIME=`date +%s`


# Run AKscore2
    $python $Akscore2 -r ${protein_pdb} -l 7e2e_C1V_gpuR_docked.pdb -s akscore2_dockc -o ./7e2e_C1Vak2.csv --batch_size 64 --num_workers 4 --device gpu



END_TIME=`date +%s`
ELAPSED_TIME=$((END_TIME-START_TIME))
echo "Total elapsed time: $ELAPSED_TIME sec"




echo $(date)
