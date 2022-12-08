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


labels = ['lambda1', 'lambda2', 'lambda_tilde', 'mass1', 'mass2', 'radius', 'mchirp', 'eos', 'radius_1p4']
samples_prior_rec2nsat = {l: [] for l in labels}
distance = 40.7

dist_eos = Uniform(eos=(1, 2000.9))
dist_m = Uniform(mass1=(1.0, 2.0), mass2=(1.0, 2.0))
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

    MRLname = 'poly5MRL/'

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


samples_prior_uniform = {l: [] for l in labels}

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

    MRLname = 'uniform/'

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

    samples_prior_uniform['mass1'].append(m1)
    samples_prior_uniform['mass2'].append(m2)
    #samples4['mchirp'].append(mchirp)
    #samples4['radius'].append(r)
    #samples4['lambda_s'].append(l)
    samples_prior_uniform['lambda1'].append(l1)
    samples_prior_uniform['lambda2'].append(l2)
    samples_prior_uniform['lambda_tilde'].append(lt)
    samples_prior_uniform['radius_1p4'].append(r_1p4)

from pycbc.inference.io import loadfile
from pycbc.conversions import lambda_tilde

samples_post_rec2nsat = {l: [] for l in labels}

fp=loadfile('/global/cscratch1/sd/bkingast/LANL_project/pycbc_runs/data/poly5_0905_ind.hdf','r+')

#parameters = fp['samples'].keys()
#samples = fp.read_samples(parameters, flatten=True)
m1=primary_mass(fp['samples']['mass1'][0:8000], fp['samples']['mass2'][0:8000])
m2=secondary_mass(fp['samples']['mass1'][0:8000], fp['samples']['mass2'][0:8000])
l1=fp['samples']['lambda1'][0:8000]
l2=fp['samples']['lambda2'][0:8000]
r1p4=fp['samples']['radius_1p4'][0:8000]

#m1=samples['mass1']
samples_post_rec2nsat['mass1'] = m1
samples_post_rec2nsat['mass2'] = m2
samples_post_rec2nsat['lambda1'] = l1
samples_post_rec2nsat['lambda2'] = l2
samples_post_rec2nsat['lambda_tilde'] = lambda_tilde(m1,m2,l1,l2)
samples_post_rec2nsat['radius_1p4'] = r1p4


samples_post_uniform = {l: [] for l in labels}

fp=loadfile('/global/cscratch1/sd/bkingast/LANL_project/pycbc_runs/data/uniform_1017_ind.hdf','r+')

#parameters = fp['samples'].keys()
#samples = fp.read_samples(parameters, flatten=True)
m1=primary_mass(fp['samples']['mass1'][0:8000], fp['samples']['mass2'][0:8000])
m2=secondary_mass(fp['samples']['mass1'][0:8000], fp['samples']['mass2'][0:8000])
l1=fp['samples']['lambda1'][0:8000]
l2=fp['samples']['lambda2'][0:8000]
r1p4=fp['samples']['radius_1p4'][0:8000]

#m1=samples['mass1']
samples_post_uniform['mass1'] = m1
samples_post_uniform['mass2'] = m2
samples_post_uniform['lambda1'] = l1
samples_post_uniform['lambda2'] = l2
samples_post_uniform['lambda_tilde'] = lambda_tilde(m1,m2,l1,l2)
samples_post_uniform['radius_1p4'] = r1p4

# plot lambda_tilde priors and posterior
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
params = ['lambda_s']
ax_labels = [r'$\tilde \Lambda$']
#ax_bounds = [(0, 5000), (0, 5000), (0, 5000), (1.6, 18)]
#for a, p, l, b in zip(ax.flatten(), params, ax_labels, ax_bounds):
#for a, p, l in zip(ax.flatten(), params, ax_labels):
#max_param_val = int(max(samples4[p]))
#num_bins = max_param_val/100 if max_param_val > 2 else 50
num_bins = 50
#ax.hist(samples_prior_unilam['lambda_tilde'], bins=40, histtype='step',
#        facecolor='None', edgecolor='red', ls='dotted', lw=2,
#            density=True, label=r"Prior: Uniform $\tilde \Lambda$")
# ax.hist(samples_post_un['lambda_tilde'], bins=num_bins, histtype='step',
#         facecolor='None', edgecolor='red', ls='solid', lw=2,
#         density=True, label=r"Posterior: Uniform $\tilde \Lambda$ prior")
ax.hist(samples_prior_rec2nsat['lambda_tilde'], bins=40, histtype='step',
        facecolor='None', edgecolor='g', ls='dotted', lw=2, density=True, label=r"Prior:poly5")
ax.hist(samples_post_rec2nsat['lambda_tilde'], bins=num_bins, histtype='step',
        facecolor='None', edgecolor='g', ls='solid', lw=2,
        density=True, label=r"Posterior:poly5 prior")
#ax.set(xlabel=l)
#         density=True, label=r"Posterior: Uniform $\tilde \Lambda$ prior")
ax.hist(samples_prior_uniform['lambda_tilde'], bins=40, histtype='step',
        facecolor='None', edgecolor='b', ls='dotted', lw=2, density=True, label=r"Prior:Uniform $R_{1.4}$")
ax.hist(samples_post_uniform['lambda_tilde'], bins=num_bins, histtype='step',
        facecolor='None', edgecolor='b', ls='solid', lw=2,
        density=True, label=r"Posterior:Uniform $R_{1.4}$ prior") 
#ax.xaxis.label.set_size(22)
#ax.yaxis.label.set_size(22)

plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel(r"$\tilde \Lambda$", fontsize=18)
plt.ylabel(r"Probability Density", fontsize=18)
plt.legend(fontsize=15)
plt.title('Comparison of cs5 and Uniform', fontsize=18, pad=18)
plt.tight_layout()
plt.savefig("../plots/poly5_lambda_prior_post.png")

# plot lambda_tilde priors and posterior
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
#ax_bounds = [(0, 5000), (0, 5000), (0, 5000), (1.6, 18)]
#for a, p, l, b in zip(ax.flatten(), params, ax_labels, ax_bounds):
#for a, p, l in zip(ax.flatten(), params, ax_labels):
#max_param_val = int(max(samples4[p]))
#num_bins = max_param_val/100 if max_param_val > 2 else 50
num_bins = 50

ax.hist(samples_prior_rec2nsat['radius_1p4'], bins=40, histtype='step',
        facecolor='None', edgecolor='g', ls='dotted', lw=2, density=True, label=r"Prior:poly5")
ax.hist(samples_post_rec2nsat['radius_1p4'], bins=num_bins, histtype='step',
        facecolor='None', edgecolor='g', ls='solid', lw=2,
        density=True, label=r"Posterior:poly5 prior")

ax.hist(samples_prior_uniform['radius_1p4'], bins=40, histtype='step',
        facecolor='None', edgecolor='b', ls='dotted', lw=2, density=True, label=r"Prior:Uniform $R_{1.4}$")
ax.hist(samples_post_uniform['radius_1p4'], bins=num_bins, histtype='step',
        facecolor='None', edgecolor='b', ls='solid', lw=2,
        density=True, label=r"Posterior:Uniform $R_{1.4}$ prior")
#ax.set(xlabel=l)
#ax.xaxis.label.set_size(22)
#ax.yaxis.label.set_size(22)

plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel(r"$R_{1.4}$ (km)", fontsize=18)
plt.ylabel(r"Probability Density", fontsize=18)
plt.legend(fontsize=15, loc="upper left")
plt.title('Comparison of lin5 and Uniform', fontsize=18, pad=18)
plt.tight_layout()
plt.savefig("../plots/poly5_r1p4_prior_post.png")
