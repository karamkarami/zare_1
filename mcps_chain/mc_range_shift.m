function [shift, rho] = mc_range_shift(pMine, pRef, cells, maxShift)
%MC_RANGE_SHIFT Range shift that lines up two mean range profiles.
%
%   [shift, rho] = mc_range_shift(pMine, pRef, cells, maxShift)
%
%   pMine, pRef : mean power along range (1 x cells)
%   cells       : cells of pRef that are used (e.g. the short-pulse part only)
%   shift       : pRef(c) best matches pMine(c + shift), |shift| <= maxShift
%   rho         : correlation of the two dB profiles at that shift

x = toDb(pMine(:));
y = toDb(pRef(:));
n = numel(x);
minOverlap = max(20, round(numel(cells) / 2));
shift = 0;
rho   = -Inf;
for s = -maxShift:maxShift
    c = cells(cells + s >= 1 & cells + s <= n);
    if numel(c) < minOverlap
        continue
    end
    r = corrOf(x(c + s), y(c));
    if r > rho
        rho   = r;
        shift = s;
    end
end
end

function d = toDb(p)
p = double(p);
d = 10*log10(max(p, max(max(p) * 1e-12, realmin)));
end

function r = corrOf(a, b)
a = a - mean(a);
b = b - mean(b);
r = sum(a .* b) / sqrt(sum(a.^2) * sum(b.^2));
end
