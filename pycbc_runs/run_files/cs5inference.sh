#!/bin/sh

# configuration files
PRIOR_CONFIG=/global/cscratch1/sd/bkingast/EOS_inference/ini_files/gw170817_cs5.ini
DATA_CONFIG=/global/cscratch1/sd/bkingast/EOS_inference/ini_files/data.ini
SAMPLER_CONFIG=/global/cscratch1/sd/bkingast/EOS_inference/ini_files/emcee_pt-gw170817.ini

OUTPUT_PATH=cs5_0829.hdf

# the following sets the number of cores to use; adjust as needed to
# your computer's capabilities

# run sampler
# Running with OMP_NUM_THREADS=1 stops lalsimulation
# from spawning multiple jobs that would otherwise be used
# by pycbc_inference and cause a reduced runtime.
OMP_NUM_THREADS=1
pycbc_inference --verbose \
    --seed 1897234 \
    --config-file ${PRIOR_CONFIG} ${DATA_CONFIG} ${SAMPLER_CONFIG} \
    --output-file ${OUTPUT_PATH} \
    --nprocesses ${SLURM_NTASKS} \
