function [main_stat, perm_stat, cfg_stats] = grp_permutations(dat, cfg_stats, verbose)
% Generate permutations to be used later by other functions in the toolbox.
%
% Recieves data which is dim1 x dim2 x subject (e.g. time x frequency, can 
% be just time x 1, time x time for temporal generalization matrices etc). 
% cfg_stats needs to have the fields grp_stat (see options below), n_perm 
% (number of permutations) and firstlevel_type (for the cluster statistic).
% 
% wilcoxon test uses MVPA-Light's implementation which should be much
% faster: https://github.com/treder/MVPA-Light.
%
% Written by Gal Vishne, Deouell Lab ~2022
% Bug reports \ requests: gal.vishne@gmail.com

if ~exist('verbose','var'); verbose = true; end
n_subj = size(dat,3);
if cfg_stats.firstlevel_type=="p"
    cfg_stats.firstlevel_type = "stat";
    switch cfg_stats.grp_stat
        case 'ttest', inv_func = @(p) tinv(p, n_subj-1); % 2nd arg is df
        case 'wilcoxon', inv_func = @(p) norminv(p);
    end
    if cfg_stats.n_sides == 1
        cfg_stats.firstlevel_thresh = inv_func(1-cfg_stats.firstlevel_thresh);
    else
        cfg_stats.firstlevel_thresh = [inv_func(cfg_stats.firstlevel_thresh/2),inv_func(1-cfg_stats.firstlevel_thresh/2)];
        % because the next steps assume pos is first!
        cfg_stats.firstlevel_thresh = fliplr(cfg_stats.firstlevel_thresh);
    end
end

perm_stat = nan([size(dat,[1 2]), cfg_stats.n_perm]);
for p = 0:cfg_stats.n_perm
    if p == 0; perm_dat = dat; else; perm_dat = dat.*sign(randn(1, 1, n_subj)); end
    switch cfg_stats.grp_stat
        case 'mean', cur_stat = mean(perm_dat, 3);
        case 'ttest', [~,~,~,sts] = ttest(perm_dat,zeros(size(perm_dat)),'Dim',3); cur_stat = sts.tstat; % zeros(size(dat)) because of a bug in Matlab input parsing
        case 'wilcoxon', cur_stat = permute(mv_stat_wilcoxon_signrank(permute(perm_dat,[3 1 2])),[2 3 1]); % MVPA light implementation that should be much faster
    end
    if p == 0; main_stat = cur_stat;
    else; perm_stat(:,:,p) = cur_stat; if verbose && mod(p,500)==0; fprintf('Done with permutation %d/%d\n', p, cfg_stats.n_perm); end
    end
end
end