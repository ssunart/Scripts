#!/bin/bash
#SBATCH --partition=rome_short    # cpu1
# #SBATCH --partition=skylake_short # cpu[2-7]
#SBATCH --ntasks=6
#SBATCH --mem=0
#SBATCH --job-name=sort_adGPU_p1
#SBATCH --output=sort_adGPU01.out
