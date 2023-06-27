#! /usr/bin/env python -v

import warnings
warnings.filterwarnings('ignore')
import numpy

from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'stixgeneral'
rcParams['mathtext.fontset'] = 'stix'

import pycbc
from pycbc import conversions, cosmology
from pycbc.distributions import Uniform, Gaussian
from pycbc.conversions import primary_mass, secondary_mass

MRL_labels = ['cs3MRL','cs5MRL','poly5MRL','lin5MRL','uniform']

labels = ['lambda1', 'lambda2', 'lambda_tilde', 'mass1', 'mass2', 'radius', 'mchirp', 'eos', 'radius_1p4']
samples_prior_rec2nsat = {l: [] for l in labels}
distance = 40.7

dist_eos = Uniform(eos=(1, 2000.9))
dist_m = Uniform(mass1=(1.0, 2.0), mass2=(1.0, 2.0))

fig, ax = plt.subplots(1,1, figsize=(8,6))
num_bins=40

for j in range(len(MRL_labels)):
    for (e,), (m1, m2) in zip(dist_eos.rvs(size=1000000), dist_m.rvs(size=1000000)):
        #if m1 < m2:
        #    continue
        m1p = primary_mass(m1, m2)
        m2s = secondary_mass(m1, m2)
        # apply mchirp
        mchirp = conversions.mchirp_from_mass1_mass2(m1p, m2s)
        if not 1.1876 < mchirp < 1.2076:
            continue

        e = int(e)

        MRLname = MRL_labels[j] + "/"

        tov_file = '/global/cscratch1/sd/bkingast/EOS_inference/EOS/LANL_Project_eos/' + MRLname + str(e) + '.dat'
        data = numpy.loadtxt(tov_file)
        radius_from_file = data[:, 0]
        mass_from_file = data[:, 1]
        lambda_from_file = data[:, 2]

        m1_src = m1p/(1.0 + pycbc.cosmology.redshift(distance))
        r1 = numpy.interp(m1_src, mass_from_file, radius_from_file)
        l1 = numpy.interp(m1_src, mass_from_file, lambda_from_file)

        m2_src = m2s/(1.0 + pycbc.cosmology.redshift(distance))
        l2 = numpy.interp(m2_src, mass_from_file, lambda_from_file)
        r2 = numpy.interp(m2_src, mass_from_file, radius_from_file)

        lt = conversions.lambda_tilde(m1_src, m2_src, l1, l2)

        r_1p4 = numpy.interp(1.4, mass_from_file, radius_from_file)

        samples_prior_rec2nsat['mass1'].append(m1)
        samples_prior_rec2nsat['mass2'].append(m2)
    #samples4['mchirp'].append(mchirp)
    #samples4['radius'].append(r)
    #samples4['lambda_s'].append(l)
        samples_prior_rec2nsat['lambda1'].append(l1)
        samples_prior_rec2nsat['lambda2'].append(l2)
        samples_prior_rec2nsat['lambda_tilde'].append(lt)
        samples_prior_rec2nsat['radius_1p4'].append(r_1p4)

#plot r1p4 priors
    ax.hist(samples_prior_rec2nsat['radius_1p4'], bins=40, histtype='step',
        facecolor='None', edgecolor='blue', ls='dotted', lw=2, density=True, label=MRL_labels[j])

plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel(r"$R_{1.4}$ (km)", fontsize=18)
plt.ylabel(r"Probability Density", fontsize=18)
plt.legend(fontsize=15)
plt.tight_layout()
plt.savefig("../plots/all_r1p4_priors.png")
