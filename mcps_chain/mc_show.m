function mc_show(X, ttl, clim)
%MC_SHOW Image of a lane in dB: range cell x row (3-D lanes: max over the bins).
%
%   mc_show(X, ttl)          automatic colour range (median .. 99.9 %; CFAR maps:
%                            all detected values)
%   mc_show(X, ttl, clim)    fixed colour range [dB]
%
%   Large lanes are thinned for the screen by keeping the maximum of each block
%   of cells / rows (at most 600 cells x 800 rows), so single detections stay visible.
%   Zeros (no detection after the CFAR) are drawn at the bottom of the range.

if ndims(X) == 3
    X = max(X, [], 3);
end
[X, xCells, yRows] = mc_thin(abs(X), 800, 600);
v = 20*log10(double(X));
if nargin < 3 || isempty(clim)
    clim = mc_clim(v);
end
v(v == -Inf) = clim(1);
imagesc(xCells, yRows, v, clim);
axis xy;
colorbar;
title(ttl, 'Interpreter', 'none');
xlabel('range cell');
ylabel('row (pulse)');
end
