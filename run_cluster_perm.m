function [cluster_masks, cluster_pvalues, cluster_stats] = run_cluster_perm(cfg_stats, main_stat, perm_stat, main_p, perm_p, verbose)
% Code for running cluster permutations according to Maris & Oostenveld, 2007
%       https://doi.org/10.1016/j.jneumeth.2007.03.024
%   (implementation of cluster mass test, when refering to 'cluster_stat'
%   this means the sum of values included in a specific cluster)
%
% Please cite: Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% The code was written as part of the analysis for this paper.
% 
% This was written originally for single subject decoding results, BUT
% it's quite general and *you can also use it for other things* as long as
% you write code to generate the permutation results outside of this
% function (and also calculate the p-values per permutation if needed).
% E.g. for the classical case of comparing means between two conditions you
% would use main_stat=t_value per time_point (and main_p accordingly) and
% the permutations would be computed by shuffling condition labels and
% running t-tests on each of these shuffles (as explained in the original
% paper). chance_level in this case whould be zero (no condition
% differences) [see below for more details on how chance_level is used].
%
% Input:
%   cfg_stats - Settings structure, see fields below.
%   main_stat - The real data result (e.g. your decoding accuracies for a
%                 single subject or t-test values for a group test)
%                   Size: n_time1 x n_time2 (n_time1\n_time2 can be == 1)
%   perm_stat - Permutation results - generate them before this function!
%                   Size: n_time1 x n_time2 x n_perm (n_perm must be 3rd! also if n_time1\n_time2 = 1)
%   main_p, perm_p - Optional, needed only if cfg_stats.firstlevel_thresh = 'p'
%                 This is single time point p_values corresponding to
%                 main_stat and perm_stat (same sizes). For the decoding
%                 case this is relevant if you used binomial test for
%                 accuracies, but it can also be used for other contexts of
%                 cluster permutations (including group level statistics)
%   verbose   - true\false (default: true) - progress updates printed or 
%                 not [it's pretty fast anyway, so currently I just print
%                 every 1000 permutations]. if you want to change this and
%                 not insert p_values use [] for those fields.
% 
% Output: 
%   cluster_masks - Cell array containing the locatins of each significant cluster
%       (size main_stat true for points in the cluster, false otherwise)
%        * Clusters are ordered from largest (largest cluster mass) to smallest
%        * Maximal number of clusters depends on cfg_stats.n_clusters
%   cluster_pvalues - corresponding array of p-values
%   cluster_stats - corresponding array of cluster mass statistics (the sum
%        of main_stat in the points belonging to the cluster)
%        * chance_level is subtracted from main_stat before this is calculated
%        * correct sign (indicates if the cluster is larger\smaller than chance)
%
%
% cfg_stats fields:
%   p_thresh          - Default: 0.05
%                         * This is for the cluster threshold \ fdr, not the cluster first level. 
%   n_sides           - Default: 1 (right side)
%                         Options: 1\2. Decides if we check only higher 
%                         from chance (right sided test) or double sided.
%   n_clusters        - Default: Inf (all clusters). Max # of clusters to return. 
%                         In the permutations only the largest cluster is taken, but we 
%                         can choose to compare more than just our largest cluster
%                         (it's just stricter, but valid, see the original paper for more details)
%   chance_level      - Default: 0. 
%                         *** Super important to adjust this if you're looking at something else!
%                         (e.g. decoding results). We need this to split to positive\negative 
%                         clusters and we assume that firstlevel_thresh was given as
%                         distance from this value if firstlevel_type == 'stat'. 
%   firstlevel_type   - No default.
%                         How should the firstlevel_thresh be treated?
%                         Options: 'p'\'stat'\'zscore'
%                           'p' - p_value -> MUST insert main_p, perm_p for this!
%                           'stat' - apply the threshold to the main_stat directly
%                           'zscore' - zscore the values according to the permutations
%                              mean & std and apply the threshold according to that.
%                              Note we assume that chance_level moves to 0 now (see below)
%   firstlevel_thresh - No default.
%                        * Note for firstlevel_type == 'stat': Insert this as DISTANCE FROM CHANCE.
%                          (e.g. if you want all auc>=0.6 to be included and you
%                          inserted chance_level = 0.5 use 0.1 as the firstlevel threshold).
%                        * For firstlevel_type == 'stat'\'zscore':
%                          If n_sides == 2 this could be 2-entry array
%                          [pos, neg] <-POSITIVE FIRST, neg really has to
%                          be negative we're taking it relative to chance_level.
%                          If only one given it's applied as [+thresh,-thresh]
%                          (in both cases it's around chance_level)
%
% Written by Gal Vishne, Deouell Lab 2019-2023
% Bug reports \ requests: gal.vishne@gmail.com

n_perm = size(perm_stat,3); % note to self: this used to come from cfg_stats, but there is really no need for that
% verify inputs \ adjust values to defaults:
if ~exist('main_p','var'); main_p = nan(size(main_stat)); end
if ~exist('perm_p','var'); perm_p = nan(size(perm_stat)); end
if ~exist('verbose','var'); verbose = true; end
if strcmp(cfg_stats.firstlevel_type,'p') && (isempty(main_p) || isempty(perm_p))
    error('No p-values input. Change firstlevel_type or add p-values.')
end
cfg_stats = cfg_stats_parser(cfg_stats);

% this is essentially a specific case of stat, so we adjust it here to simplify things
if strcmp(cfg_stats.firstlevel_type,'zscore')
    cfg_stats.firstlevel_type = 'stat'; cfg_stats.chance_level = 0; % this is the assumption after z-scoring (might be best to use chance_level instead of the perm_mean in the future)
    perm_mean = mean(perm_stat,3); perm_std = std(perm_stat,[],3);
    main_stat = (main_stat - perm_mean)./perm_std;
    perm_stat = (perm_stat - perm_mean)./perm_std;
end

[cStats, cMasks] = find_clusters(main_stat, main_p, cfg_stats);
% Sort clusters largest to smallest
[~,c_sort_idx] = sort(cStats,'descend');
cStats = cStats(c_sort_idx);
cStats = cStats(1:min(cfg_stats.n_clusters, length(cStats))); % take out only the clusters we asked for

permutation_distribution = nan(n_perm,1); % get the cluster statistic from the largest cluster of each permutation
for p = 1:n_perm
    cluster_stats_perm = find_clusters(perm_stat(:,:,p), perm_p(:,:,p), cfg_stats);
    if isempty(cluster_stats_perm)
        permutation_distribution(p) = 0;
    else
        permutation_distribution(p) = max(cluster_stats_perm);
    end
    if mod(p,1000)==0 && verbose; fprintf('Done with cluster stats for permutation %d/%d\n',p,n_perm); end
end

cluster_masks = cell(size(cStats)); % turn it from one array with numbered clusters to a cell array with each cluster in a separate cell
cluster_pvalues = nan(size(cStats));
cluster_stats = nan(size(cStats));
for c_idx = 1:length(cStats)
    cluster_pvalues(c_idx) = (sum(permutation_distribution >= cStats(c_idx))+1) / (n_perm+1); % In find_clusters we turned all cStats and all permutation_distribution to positive clusters, so this is fine also if you are doing a 2-sided test
    cluster_masks{c_idx}   = (cMasks == c_sort_idx(c_idx));
    cluster_stats(c_idx)   = cStats(c_idx)*sign(sum(main_stat(cluster_masks{c_idx})-cfg_stats.chance_level,'all')); % Fix back the sign of the cluster_stats (relative to chance_level)
end

% Return only significant clusters
signif_idx = cluster_pvalues<=cfg_stats.p_thresh;
if any(signif_idx)
    cluster_masks = cluster_masks(signif_idx);cluster_stats = cluster_stats(signif_idx);cluster_pvalues = cluster_pvalues(signif_idx);
else
    cluster_masks = {zeros(size(main_stat))};cluster_pvalues=[];cluster_stats = [];
end
end

function [cStats, cMasks] = find_clusters(mat_stat, mat_p, cfg_stats)
% Input:
% - mat_stat - 1d\2d array
% - mat_p only used if firstlevel_type=='p'
% - cfg_stats needed fields: firstlevel_type, firstlevel_thresh, n_sides, chance_level
% REMEMBER: if firstlevel_type=='stat'\'zscore' and n_sides==2 and you want to
%   insert 2 thresholds write them as [POS,NEG] (sorry about that, 
%   backwards compatibility for me...). If only one value is given it will
%   be applied as chance+[THRESH,-THRESH].
% Output:
% - cStats - array with the cluster stats -> turns everything to positive
%     (we fix that later in the main function)
% - cMasks - array the size of mat_stat, with numbers 0:n_clusters. 0 means
%     the point does not belong to any cluster, 1:n_clusters mark which
%     cluster the point belongs to (numbers matched to indices in cStats)
chance_level = cfg_stats.chance_level; firstlevel_thresh = cfg_stats.firstlevel_thresh;
switch cfg_stats.firstlevel_type
    case 'p'
        sig_all = mat_p <= firstlevel_thresh;
        sig_pos = sig_all & mat_stat > chance_level;
        sig_neg = sig_all & mat_stat < chance_level;
    case 'stat'
        sig_pos = mat_stat >= firstlevel_thresh(1) + chance_level;
        if cfg_stats.n_sides == 2
            if numel(firstlevel_thresh) == 2
                sig_neg = mat_stat <= firstlevel_thresh(2) + chance_level;
            else
                sig_neg = mat_stat <= -firstlevel_thresh + chance_level;
            end
        end
    otherwise
        error("Only recognizes firstlevel_thresh = 'p'\'stat' at this point, retrace your steps")
end
cMasks = zeros(size(mat_stat));
% Find the above-threshold clusters: (positive clusters)
%   Uses the bwconncomp function in Matlab which finds clusters in binary
%   images. For more details see the Matlab documentation of the function.
CC1 = bwconncomp(sig_pos, 4); % 4 means no diagonal connections
cStats = nan(1, CC1.NumObjects);
for clust=1:CC1.NumObjects
    cMasks(CC1.PixelIdxList{clust}) = clust;
    cStats(clust) = sum(mat_stat(CC1.PixelIdxList{clust})-chance_level,'all');
end
if cfg_stats.n_sides == 2 % If relevant, also look for negative clusters
    CC2 = bwconncomp(sig_neg,4);
    cStat2 = nan(1, CC2.NumObjects);
    for clust=1:CC2.NumObjects
        cMasks(CC2.PixelIdxList{clust}) = CC1.NumObjects+clust;
        cStat2(clust) = sum(chance_level-mat_stat(CC2.PixelIdxList{clust}),'all'); % ! Make all cluster stats positive (important to determine order & who is maximal later)
    end
    cStats = [cStats,cStat2];
end
end

function cfg_stats = cfg_stats_parser(cfg_stats)
% Make sure all fields are there & correct. Add default values for missing things.
if ~isfield(cfg_stats,'p_thresh')
    warning('Setting p_thresh to 0.05'); cfg_stats.p_thresh = 0.05;
end
if ~isfield(cfg_stats,'n_sides')
    warning('Setting n_sides to 1 (larger than, right sided only)'); cfg_stats.n_sides = 1;
end
if ~isfield(cfg_stats,'n_clusters')
    warning('Setting n_clusters to Inf (all ordered from largest to smallest)'); cfg_stats.n_clusters = Inf;
end
if ~isfield(cfg_stats,'chance_level')
    warning("Setting chance_level to 0. If your chance level is different it's important to change it!"); cfg_stats.chance_level = 0;
end
if ~isfield(cfg_stats,'firstlevel_type') 
    error('You must set field: firstlevel_type');
end
if ~isfield(cfg_stats,'firstlevel_thresh')
    error('You must set field: firstlevel_thresh');
end
end