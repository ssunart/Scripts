##AMBER 기반의 기본적인 MD 워크플로우 
#SLURM 작업 스케줄러를 사용하여 AMBER 프로그램으로 분자 동역학(MD) 시뮬레이션을 실행하는 스크립트

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 28
##SBATCH --cpus-per-task=5
##SBATCH -p gpu_short
#SBATCH -p skylake_short
##SBATCH --gres=gpu:1
#SBATCH --job-name=MD_NNK
#SBATCH --output=MD_NNK.out
#sbatch dock.sh < Input your command console

echo "Start Time = `date`"

module load compiler/gcc-7.5.0 cuda/11.3 mpi/openmpi-4.0.2-gcc-7
source /home/juyong/apps/amber20_3090/amber.sh


amber=pmemd.cuda
init=step3_input
mini_prefix=step4.0_minimization
equi_prefix=step4.1_equilibration
prod_prefix=step5_production
prod_step=step5

#${amber} -O -i ${mini_prefix}.mdin -p ${init}.parm7 -c ${init}.rst7 -o ${mini_prefix}.mdout -r ${mini_prefix}.rst7 -inf ${mini_prefix}.mdinfo -ref ${init}.rst7
mpirun -np 10 pmemd.MPI -O -i ${mini_prefix}.mdin -p ${init}.parm7 -c ${init}.rst7 -o ${mini_prefix}.mdout -r ${mini_prefix}.rst7 -inf ${mini_prefix}.mdinfo -ref ${init}.rst7

exit
${amber} -O -i ${equi_prefix}.mdin -p ${init}.parm7 -c ${mini_prefix}.rst7 -o ${equi_prefix}.mdout -r ${equi_prefix}.rst7 -inf ${equi_prefix}.mdinfo -ref ${init}.rst7 -x ${equi_prefix}.nc

cnt=1
cntmax=1
while (( ${cnt} <= ${cntmax} )); do
    let pcnt=${cnt}-1
    istep=${prod_step}_${cnt}
    pstep=${prod_step}_${pcnt}

    if [[ ${cnt} == 1 ]]; then
	pstep=${equi_prefix}
    fi

    ${amber} -O -i ${prod_prefix}.mdin -p ${init}.parm7 -c ${pstep}.rst7 -o ${istep}.mdout -r ${istep}.rst7 -inf ${istep}.mdinfo -x ${istep}.nc
    let cnt=${cnt}+1
done
