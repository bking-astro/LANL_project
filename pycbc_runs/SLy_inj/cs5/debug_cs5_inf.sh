#!/bin/bash
#SBATCH --qos=debug
#SBATCH --constraint=cpu
#SBATCH --time=00:30:00                                                        
#SBATCH --nodes=4 
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=king0690@umn.edu
#SBATCH --mail-type=ALL                                                        
#SBATCH --output=slurm-%J.out
#SBATCH --error=slurm-%J.err

# load all the modules needed to run mpi4py
module load python
module load evp-patch
export MPI4PY_RC_RECV_MPROBE=0
export CXI_FORK_SAFE=1
export CXI_FORK_SAFE_HP=1
module load fast-mkl-amd

cd /pscratch/sd/b/bkingast/inj_inference/SLy_inj/cs5
mkdir -p logs_mpi
conda activate venv_pycbc

#export HDF5_USE_FILE_LOCKING=FALSE
export MPI4PY_RC_RECV_MPROBE=0

echo ${SLURM_NTASKS}

echo "Start time:"
date

srun -n ${SLURM_NTASKS} cs5_inf.sh 1>logs_mpi/cs5_SLy_debug.out 2>logs_mpi/cs5_SLy_debug.err

echo "End time:"
date


