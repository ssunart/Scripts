#!/bin/bash
#SBATCH --job-name=VS_test
#SBATCH --ntasks=4
#SBATCH --partition=g3090_veryshort
#SBATCH --output=./%j.log
#SBATCH --gres=gpu:1

echo $(date)

module load cuda/11.8

python=/appl/anaconda3/envs/VS_env/bin/python  # python for Akscore2 & RTMScore

Akscore2=/appl/git-repo/AKScore2/run.py     # running Akscore2
RTMScore=/appl/git-repo/RTMScore/rtmscore.py         # running RTMScore
 
mol2pdb=/appl/git-repo/AKScore2/convert_mol2pdb.py   # Making input ligand PDB for Akscore

protein_pdb='./7S15_protein.pdb'     # INPUT protein pdb for Akscore2 & RTMScore
ligand_pdb='./7S15_docked.pdb'     # INPUT reference ligand file (sdf, mol2) for RTMScore / (if you don't have pocket)

# Run AKscore2
#$python $mol2pdb -l ${ligand_mol2} -o ${ligand_pdb}
$python $Akscore2 -r ${protein_pdb} -l ${ligand_pdb} -s akscore2_dockc -o ./result_pdb.csv --batch_size 64 --num_workers 4 --device gpu

# Run RTMScore with gen pocket
#$python $RTMScore -p ${protein_pdb} -l ${ligand_sdf} -gen_pocket -c 10.0 -rl './ligand.mol2' -o 'with_pocket'

# Run RTMScore with pocket
#$python $RTMScore -p ${pocket_pdb} -l ${ligand_mol2} -o 'without_pocket'

echo $(date)
