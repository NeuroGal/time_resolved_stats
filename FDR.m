function [ ind, thres ] = FDR( p_list, alpha, corrected, num_tests )
% Computes the False Discovery Rate according to Benjamini and Hochberg (1995). 
% 
% Inputs: 
% p_list - list of p values
% alpha - the desired alpha threshold. Default: 0.05
% corrected - set to true if correction for dependencies is to be applied, according to Benjamini
% and Yekutieli (2001) (this is probably not the common case). 
% num_tests - the number of comparisons - in case the given p-values are
% only a subset of the total number of tests. 
%
% outputs:
% ind - the indexes of significant p-values within p_list
% thres - the p-value which served as the actual threshold in this test. 
% 
% Written by Edden Gerber, lab of Leon Y. Deouell, 2012
% Please send bug reports and requsts to edden.gerber@gmail.com
%
n_vals = length(p_list);

if size(p_list,1) == 1
    p_list = p_list';
end

if nargin < 2
    alpha = 0.05;
end

if nargin < 3
    corrected = false;
end
if nargin < 4
    num_tests = n_vals; 
end

p_sorted = sort(p_list,'descend');

if corrected
    comp = (num_tests:-1:1)/num_tests * alpha / sum((1:num_tests)/num_tests);
else
    comp = (num_tests:-1:1)/num_tests * alpha;
end


comp = comp((end-n_vals+1):end)';

i = find(p_sorted <= comp,1,'first');

if isempty(i)
    thres = 0;
else
    thres = p_sorted(i);
end

ind = find(p_list<=thres);

end

