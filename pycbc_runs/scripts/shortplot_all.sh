#!/bin/sh
pycbc_inference_plot_posterior --verbose \
        --input-file ../data/poly5_0905_ind.hdf ../data/uniform_1017_ind.hdf \
        --output-file shortpoly5uniform_wprior.png \
        --plot-contours \
        --plot-marginal \
	--plot-prior /global/cscratch1/sd/bkingast/EOS_inference/ini_files/gw170817_poly5.ini /global/cscratch1/sd/bkingast/EOS_inference/ini_files/gw170817_uniform.ini \
	--parameters \
"mchirp_from_mass1_mass2(mass1, mass2)/(1+redshift(40.7)):$\mathcal{M}^{src} (\mathrm{M}_{\odot})$" \
"lambda_tilde(mass1, mass2, lambda1, lambda2):\$\tilde{\Lambda}$" \
"radius_1p4:\$R_{1.4} (km)$" \
	--input-file-labels "poly5" "uniform" \
        --iteration -1
