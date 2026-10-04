function v = mc_prctile(x, p)
%MC_PRCTILE Percentiles of the finite values of x (no toolbox), on at most 1e6 values.
%
%   v = mc_prctile(x, [50 99.9])

x = double(x(isfinite(x)));
if isempty(x)
    v = zeros(size(p));
    return
end
if numel(x) > 1e6
    x = x(round(linspace(1, numel(x), 1e6)));
end
x = sort(x);
v = x(min(max(round(p / 100 * numel(x)), 1), numel(x)));
v = reshape(v, size(p));
end
