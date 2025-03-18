# time_resolved_stats
Matlab toolbox with functions for performing statistical quantification and correction for multiple comparisons using permutation-based methods (intended as correction across time-points, but you can use it also for time x frequency data and more).

This was written as part of the analysis for the following manuscript (so please cite if you use it):
Vishne et al., Cell Reports 2023, 'Distinct Ventral Stream and Prefrontal Cortex Representational Dynamics during Sustained Conscious Visual Perception' (biorxiv DOI, to be updated when formally published): https://doi.org/10.1016/j.celrep.2023.112752.

Code from this repository is used by my repositories for single-subject decoding (**iEEG_decoding_minitoolbox**, https://github.com/NeuroGal/iEEG_decoding_minitoolbox) and representational similarity analyses (**Gals_RSA_toolbox**, https://github.com/NeuroGal/Gals_RSA_toolbox).


Implemented options:
- ‘Perm_pvals.m’ - get p-values from point-by-point permutations. The p-values here are not corrected for multiple comparisons, this can be done using ‘FDR.m’ written by Edden Gerber (https://github.com/edden-gerber/time_series_analysis_and_statistics/blob/master/testing/FDR.m). This implements FDR correction according to Benjamini and Hochberg (1995). A version of this function which was on our lab server and definitely worked with this code is also enclosed.
- ‘Max_stat_correction.m’ - max-statistic correction for multiple comparisons (controlling the FWER). See Nichols & Holmes, 2002; https://doi.org/10.1002/hbm.1058 for more details.
- ‘Run_cluster_perm.m’ - cluster-based permutations correction for multiple comparisons (FWER), see Maris & Oostenveld, 2007; https://doi.org/10.1016/j.jneumeth.2007.03.024 for more details. An early version of this function is based on Edden Gerber’s code (see also https://edden-gerber.github.io/eeg_functions/).
All of these functions require you to create the permutation data outside of the function. For decoding\RSA this is done in the toolboxes mentioned above, for group-level permutations you can use the function ‘grp_permutations.m’ (for the options included see the text in the function).


Gal Vishne, June 2023

Twitter: @neuro_gal
