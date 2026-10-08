function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method)
if nargin < 2 || isempty(q), q = 0.05; end
if nargin < 3 || isempty(method), method = 'pdep'; end

p = pvals(:);
m = length(p);
[p_sorted, sort_ids] = sort(p);
orig_ids = 1:length(p);
orig_ids(sort_ids) = orig_ids;

if strcmpi(method, 'pdep')
    adj_p = m * p_sorted ./ (1:m)';
else
    adj_p = m * p_sorted ./ (1:m)' / sum(1./(1:m));
end

adj_p = cummin(adj_p(end:-1:1));
adj_p = adj_p(end:-1:1);
adj_p(adj_p > 1) = 1;
adj_p = adj_p(sort_ids);

h = adj_p <= q;
crit_p = max(p(h)); if isempty(crit_p), crit_p = 0; end
adj_ci_cvrg = 1 - crit_p;
end