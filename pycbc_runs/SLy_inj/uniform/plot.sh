#!/bin/sh
pycbc_inference_plot_posterior --verbose \
        --input-file uniform_SLy_1p41p4_ind.hdf \
        --output-file posterior_uniform_SLy_1p41p4_run3.png \
        --plot-scatter \
        --plot-contours \
        --plot-marginal \
	--parameters \
"mchirp_from_mass1_mass2(mass1, mass2)/(1+redshift(40.7)):$\mathcal{M}^{src} (\mathrm{M}_{\odot})$" \
"primary_mass(mass1,mass2)/(1+redshift(40.7)):\$m_1^{src} (\mathrm{M}_{\odot})$" \
"secondary_mass(mass1,mass2)/(1+redshift(40.7)):\$m_2^{src} (\mathrm{M}_{\odot})$" \
"primary_mass(mass1, mass2)/secondary_mass(mass1, mass2):\$q$" \
"lambda_tilde(mass1, mass2, lambda1, lambda2):\$\tilde{\Lambda}$" \
"radius1:\$R_1 (km)$" \
"radius2:\$R_2 (km)$" \
"radius_1p4:\$R_{1.4} (km)$" \
        --z-arg snr \
        --iteration -1
