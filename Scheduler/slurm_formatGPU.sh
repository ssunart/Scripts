#!/bin/bash
#SBATCH --partition=g3090_short    # gpu2
# #SBATCH --partition=a5000_short # gpu3[2-7]
# #SBATCH --partition=icelake  # cpu[8-11]
#SBATCH --ntasks=2
#SBATCH --mem=0
#SBATCH --gres=gpu:1
#SBATCH --job-name=whole_gpu7_01
#SBATCH --output=whole_adGPU_p7_01.out

