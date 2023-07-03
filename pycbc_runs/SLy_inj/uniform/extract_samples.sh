#!/bin/bash -v

path_to_original_hdf_file=uniform_SLy_1p41p4.hdf.checkpoint
path_to_independent_samples_file=uniform_SLy_1p41p4_ind.hdf

pycbc_inference_extract_samples --input-file $path_to_original_hdf_file --output-file $path_to_independent_samples_file --verbose --force

