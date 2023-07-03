#!/bin/sh

# configuration files
PRIOR_CONFIG=/pscratch/sd/b/bkingast/inj_inference/SLy_inj/cs5/gw170817_cs5.ini
DATA_CONFIG=/pscratch/sd/b/bkingast/inj_inference/SLy_inj/data_SLy_1p41p4.ini
SAMPLER_CONFIG=/pscratch/sd/b/bkingast/EOS_inference/ini_files/emcee_debug.ini

OUTPUT_PATH=cs5_SLy_1p41p4_test.hdf

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
    --nprocesses 160 \
    --use-mpi \
    --processing-scheme cpu \
