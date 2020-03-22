If you use our work in your research, please cite as:

P. Irofti, F. Stoican, and V. Puig, Fault Handling in Large Water Networks with Online Dictionary Learning, arXiv preprint arXiv:2003.08483, 2020
```
@book{ISP20_ddnet-online,
author    = {Irofti, P. and Stoican, F., and Puig, V.},
title     = {Fault Handling in Large Water Networks with Online Dictionary
Learning},
year = {2020},
eprint = {2003.08483},
archiveprefix = {arXiv}
}
```

Requirements:


o the [dictionary learning toolbox](https://github.com/pirofti/dl-box)

o the [YALMIP](https://yalmip.github.io/) solver

o for generating the data you will need the
[EPANET-Matlab-Toolkit](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit).

Main scripts:

o get_residuals_random_profile.m, it generates and saves (as a mat file) the residual data for the given network

o run_fdi_toddler.m, runs the [TODDLeR](https://github.com/abaltoiu/malid) algorithm over the pre-computed residuals, for various numbers of placed sensors; for each sensor number it saves a different mat files containing the result (success rates, dictionaries, etc)

o plot_article_figures.m, with the exception of the networks, it generates all the plots appearing in the article (some differences appear due to changes in the mat files used; nothing major)


Necessary script runs in order to obtain the figures shown in the related article


1. run get_residuals_random_profile.m with the variables containing the 'hanoi' word uncommented; the result is the residues_hanoi.mat file, saved in the ./data subfolder

2. run get_residuals_random_profile.m with the variables containing the 'hanoi_nominal' word uncommented (lines 121-122); the result is the residues_hanoi_nominal.mat file, saved in the ./data subfolder

3. run get_residuals_random_profile.m with the variables containing the 'generic_2' word uncommented; the result is the residues_generic_2.mat file, saved in the ./data subfolder

4. run run_fdi_toddler.m with the variables containing the 'generic 2' word uncommented; the result are mat files in the format data_dl_generic_2_s.mat where 's' is the replaced by the sensor number (by default, a number from 5 to 50), saved in the '/data/absolute residual generic 2/' subfolder

5. change fault_at_node_level=0 and run the previous step; the result are mat files in the format data_dl_generic_2_s_comm.mat where 's' is the replaced by the sensor number (by default, a number from 5 to 50), saved in the '/data/absolute residual generic 2 comm/' subfolder; the file ./data/generic_2-communities.mat has to exist (obtained through the python communities script)
