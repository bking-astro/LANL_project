#! /usr/bin/env python -v

import argparse
from pycbc.inference.io import loadfile
from pycbc.inference import option_utils
import numpy as np
from pycbc.conversions import primary_mass, secondary_mass
import h5py
from pycbc.cosmology import redshift

parser = argparse.ArgumentParser()
parser.add_argument("--input-file", type=str,
                    help="Inference hdf file to input.")
parser.add_argument("--verbose", action="store_true", default=False)

# add results group options
#option_utils.add_inference_results_option_group(parser)

args = parser.parse_args()

input_file = args.input_file
fp=h5py.File(args.input_file,'r+')

samples_new = {}
samples_new['mass1'] = []
samples_new['mass2'] = []
samples_new['lambda1'] = []
samples_new['lambda2'] = []
samples_new['radius1'] = []
samples_new['radius2'] = []
samples_new['radius_1p4'] = []

#for va in fp.variable_params:
#    samples[va] = fp.read_samples(va)
#    print("len(samples[va])", len(samples[va]))
#    samples_new[va] = []

for i in range(len(fp["samples"][fp.attrs['variable_params'][0]])):
    m1 = primary_mass(fp["samples/mass1"][i], fp["samples/mass2"][i])
    m2 = secondary_mass(fp["samples/mass1"][i], fp["samples/mass2"][i])
    m1_src = m1/(1.0 + redshift(40.7))
    idx1 = int(fp["samples/eos"][i])
    tov_data = np.loadtxt('/global/cscratch1/sd/bkingast/EOS_inference/EOS/LANL_Project_eos/cs3MRL/' + str(idx1) + '.dat')
    mass_data = tov_data[:,1]
    lambda_data = tov_data[:,2]
    radius_data = tov_data[:,0]
    lambdav1 = np.interp(m1_src, mass_data, lambda_data)
    r1 = np.interp(m1_src, mass_data, radius_data)
    samples_new['mass1'].append(m1)
    samples_new['lambda1'].append(lambdav1)
    samples_new['radius1'].append(r1)
    
    m2_src = m2/(1.0 + redshift(40.7))
    lambdav2 = np.interp(m2_src, mass_data, lambda_data)
    r2 = np.interp(m2_src, mass_data, radius_data)
    samples_new['mass2'].append(m2)
    samples_new['lambda2'].append(lambdav2)
    samples_new['radius2'].append(r2)

    r_1p4 = np.interp(1.4, mass_data, radius_data)
    samples_new['radius_1p4'].append(r_1p4)
    
 
#fp["samples/mass1"] = np.array(samples_new['mass1'])
#fp["samples/mass2"] = np.array(samples_new['mass2'])
if not 'lambda1' in fp['samples']:
    fp.create_dataset("samples/lambda1", data=np.array(samples_new['lambda1']))
if not 'lambda2' in fp['samples']:
    fp.create_dataset("samples/lambda2", data=np.array(samples_new['lambda2']))
if not 'radius1' in fp['samples']:
    fp.create_dataset("samples/radius1", data=np.array(samples_new['radius1']))
if not 'radius2' in fp['samples']:
    fp.create_dataset("samples/radius2", data=np.array(samples_new['radius2']))
if not 'radius_1p4' in fp['samples']:
    fp.create_dataset("samples/radius_1p4", data=np.array(samples_new['radius_1p4']))
fp.close()
print("Done")
