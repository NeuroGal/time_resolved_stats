function p_values = perm_pvals(main_stat, perm_stat, cfg_stats)
% Helper function to compute permutation p-values. 
% Inputs:
%   main_stat - The real data (e.g. your decoding accuracies for a
%               single subject or t-test values for a group test)
%               Size: n_time1 x n_time2 (n_time1\n_time2 can be == 1)
%   perm_stat - Permutation results - generate them before this function!
%               Size: n_time1 x n_time2 x n_perm (n_perm must be 3rd! also
%               if n_time1\n_time2 = 1)
%   cfg_stats - Settings structure, uses the following inputs:
%                > n_sides - 1/2: default = 1
%                   1-sided test is only the positive tail! (larger than)
%                   2-sided test checks both sides.
%                > chance_level - NO DEFAULT (mandatory for n_sides==2)
%               
% Output: point-by-point permutation values (no correction for multiple
%         comparisons)
%
% Written by Gal Vishne, Deouell Lab ~2021
% Bug reports \ requests: gal.vishne@gmail.com
% Please cite: Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% The code was written as part of the analysis for this paper.

if ~isfield(cfg_stats,'n_sides')
    warning('Setting n_sides to 1 (larger than, right sided only)'); cfg_stats.n_sides = 1;
end
if cfg_stats.n_sides == 2 && ~isfield(cfg_stats,'chance_level')
    error('You must set a value for chance_level')
end
n_perm = size(perm_stat,3);
if cfg_stats.n_sides == 1
    p_values = (sum(perm_stat>=main_stat,3)+1)/(n_perm+1);
elseif cfg_stats.n_sides == 2
    p_values = (sum(abs(perm_stat-cfg_stats.chance_level)>=abs(main_stat-cfg_stats.chance_level),3)+1)/(n_perm+1);
end
end