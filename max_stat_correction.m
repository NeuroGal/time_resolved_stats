function [max_mask, max_pvals, max_thresh] = max_stat_correction(main_stat, perm_stat, cfg_stats)
% Helper function to compute max-statistic permutations control for
% multiple comparisons - see for example: 
%   Nichols & Holmes (2002). Nonparametric permutation tests for functional
%   neuroimaging: a primer with examples. Hum. Brain Mapp. 10.1002/hbm.1058.
%
% Inputs:
%   main_stat - The real data (e.g. your decoding accuracies for a
%               single subject or t-test values for a group test)
%               Size: ntime1 x ntime2 (ntime1\ntime2 can be == 1)
%   perm_stat - Permutation results - generate them before this function!
%               Size: ntime1 x ntime2 x n_perm (n_perm must be 3rd! also
%               if ntime1\ntime2 == 1)
%       ---> i.e. we assume the permutations are in the 3rd dimension and
%            you want to correct for dimensions 1-2 (if there are other
%            dimensions afterwards they are each corrected separately)
%   cfg_stats - Settings structure, uses the following inputs:
%                > n_sides - 1/2: default = 1
%                   1-sided test is only the positive tail! (larger than)
%                   2-sided test checks both sides.
%                > p_thresh - default = 0.05
%                > chance_level - NO DEFAULT (mandatory for n_sides==2)
%               
% Output: 
%   max_mask   - true\false of significant points (ntime1 x ntime2)
%   max_pvals  - proportion of permutations each point is larger than
%                (corrected p-values, ntime1 x ntime2)
%   max_thresh - single value, the maximum correction threshold
%                * if 2-sided this is relative to chance-level (!)
%
% Written by Gal Vishne, Deouell Lab ~2022
% Bug reports \ requests: gal.vishne@gmail.com
% Please cite: Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% The code was written as part of the analysis for this paper.

dim_crct = [1 2]; dim_perm = 3;
n_perm = size(perm_stat,dim_perm);
n_sides = cfg_stats.n_sides;
if n_sides == 2
    chance_level = cfg_stats.chance_level;
    main_stat = abs(main_stat-chance_level);
    perm_stat = abs(perm_stat-chance_level);
end
max_per_perm = max(perm_stat,[],dim_crct);
max_thresh = prctile(max_per_perm, 100*(1-cfg_stats.p_thresh), dim_perm);
max_mask = max_thresh <= main_stat;
max_pvals = (sum(main_stat < max_per_perm,dim_perm)+1)/(n_perm+1);
end