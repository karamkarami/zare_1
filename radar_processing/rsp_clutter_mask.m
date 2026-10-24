function [W, azW, rW] = rsp_clutter_mask(L, azDeg, rangeM)
%RSP_CLUTTER_MASK Relative mean power of one clutter layer on an azimuth x range grid.
%
%   [W, azW, rW] = rsp_clutter_mask(L, azDeg, rangeM)
%
%   L      : one layer (rsp_clutter_types)
%   azDeg  : azimuths [deg] (column), rangeM : ranges [m] (row)
%   W      : numel(azDeg) x numel(rangeM), 0..1: soft sector edges (sin^2
%            ramps of L.edgeDeg / L.edgeM inside the region) times the slope
%            profile (0 dB at the far edge, -L.slopeDb at the near edge)
%   azW, rW: the azimuth and range factors (without the slope)
%
%   No ramp at a near edge of 0 m (the clutter starts at the radar) or for a
%   full circle.

azDeg  = azDeg(:);
rangeM = rangeM(:).';

% azimuth
a1 = L.azDeg(1);
w  = L.azDeg(2) - L.azDeg(1);
if w >= 360 || w <= -360
    azW = ones(size(azDeg));
else
    w = mod(w, 360);
    u = mod(azDeg - a1, 360);
    azW = ramp(min(u, w - u), L.edgeDeg) .* (u <= w);
end

% range
r1 = L.rangeM(1);  r2 = L.rangeM(2);
in = rangeM >= r1 & rangeM <= r2;
rW = ramp(r2 - rangeM, L.edgeM) .* in;
if r1 > 0
    rW = rW .* ramp(rangeM - r1, L.edgeM);
end
sl = 10.^(L.slopeDb * (min(max((rangeM - r1) / max(r2 - r1, eps), 0), 1) - 1) / 10);

W = azW * (rW .* sl);
end

function y = ramp(t, edge)
if edge <= 0
    y = double(t >= 0);
else
    y = sin(pi/2 * min(max(t / edge, 0), 1)).^2;
end
end
