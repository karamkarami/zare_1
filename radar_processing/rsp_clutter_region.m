function [id, names] = rsp_clutter_region(L, azDeg, rangeM, marginDeg, marginM)
%RSP_CLUTTER_REGION Which clutter region each point belongs to.
%
%   [id, names] = rsp_clutter_region(L, azDeg, rangeM)
%   [id, names] = rsp_clutter_region(L, azDeg, rangeM, marginDeg, marginM)
%
%   L       : clutter layers (rsp_clutter_types), e.g. rsp_clutter_types(S.clutter)
%   azDeg, rangeM : points (same size), [deg] and [m]
%   marginDeg, marginM : the regions are widened by this much (the beam
%             spreads clutter over about one beamwidth); default 0
%   id      : 0 = clear (noise only), k = layer k; where layers overlap the
%             strongest one; behind a mountain (in its shadow) 'shadow'
%   names   : region names, names{id + 1}: {'clear', L.name ..., 'shadow'}

if nargin < 4, marginDeg = 0; end
if nargin < 5, marginM = 0; end
sz = size(azDeg);
az = azDeg(:);
r  = rangeM(:);
id = zeros(numel(az), 1);
best = -Inf(numel(az), 1);
for k = 1:numel(L)
    w  = L(k).azDeg(2) - L(k).azDeg(1);
    if abs(w) >= 360
        inA = true(size(az));
    else
        u   = mod(az - L(k).azDeg(1) + marginDeg, 360);
        inA = u <= mod(w, 360) + 2*marginDeg;
    end
    inR = r >= L(k).rangeM(1) - marginM & r <= L(k).rangeM(2) + marginM;
    sel = inA & inR & L(k).powerDb > best;
    id(sel)   = k;
    best(sel) = L(k).powerDb;
end
names = [{'clear'}, {L.name}, {'shadow'}];
% shadow: points behind a mountain that are not themselves in a region
for k = 1:numel(L)
    if L(k).shadow
        u   = mod(az - L(k).azDeg(1), 360);
        inA = u <= mod(L(k).azDeg(2) - L(k).azDeg(1), 360);
        sel = id == 0 & inA & r > L(k).rangeM(2) + marginM;
        id(sel) = numel(L) + 1;
    end
end
id = reshape(id, sz);
end
