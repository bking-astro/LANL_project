#!/bin/sh
pycbc_inference_plot_posterior --verbose \
        --input-file ../data/cs5_0829_ind.hdf \
        --output-file shortcs5_wprior.png \
        --plot-contours \
        --plot-marginal \
	--plot-prior /global/cscratch1/sd/bkingast/EOS_inference/ini_files/gw170817_cs5.ini \
	--parameters \
"mchirp_from_mass1_mass2(mass1, mass2)/(1+redshift(40.7)):$\mathcal{M}^{src} (\mathrm{M}_{\odot})$" \
"lambda_tilde(mass1, mass2, lambda1, lambda2):\$\tilde{\Lambda}$" \
"radius_1p4:\$R_{1.4} (km)$" \
        --z-arg snr \
	--input-file-labels "cs3" "cs5" \
        --iteration -1
