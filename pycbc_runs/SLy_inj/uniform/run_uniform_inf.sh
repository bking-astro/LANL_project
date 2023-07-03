#!/bin/bash
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --time=12:00:00                                                        
#SBATCH --nodes=4 
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=king0690@umn.edu
#SBATCH --mail-type=ALL                                                        
#SBATCH --output=slurm-%J.out
#SBATCH --error=slurm-%J.err

module load python
module load evp-patch
export MPI4PY_RC_RECV_MPROBE=0
export CXI_FORK_SAFE=1
export CXI_FORK_SAFE_HP=1
module load fast-mkl-amd

cd /pscratch/sd/b/bkingast/inj_inference/SLy_inj/uniform
mkdir -p logs_mpi

conda activate venv_pycbc

export OMP_NUM_THREADS=32

echo ${SLURM_NTASKS}

echo "Start time:"
date

srun -n ${SLURM_NTASKS} uniform_inf.sh 1>logs_mpi/uniform_SLy.out 2>logs_mpi/uniform_SLy.err

echo "End time:"
date


