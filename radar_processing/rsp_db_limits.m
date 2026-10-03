function lim = rsp_db_limits(xDb, low, high)
%RSP_DB_LIMITS Colour limits from percentiles of a dB image (no toolbox needed).
%
%   lim = rsp_db_limits(xDb, low, high)     defaults: low = 50, high = 99.95

if nargin < 2, low  = 50;    end
if nargin < 3, high = 99.95; end
v = xDb(isfinite(xDb));
if isempty(v)
    lim = [0 1];
    return
end
if numel(v) > 2e6
    v = v(round(linspace(1, numel(v), 2e6)));       % thin very large images
end
v   = sort(v(:));
lim = double(v(max(1, round([low high]/100 * numel(v)))))';
if lim(2) <= lim(1)
    lim = lim(1) + [0 1];
end
end
