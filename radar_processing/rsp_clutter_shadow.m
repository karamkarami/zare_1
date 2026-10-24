function F = rsp_clutter_shadow(L, azDeg, rangeM, except)
%RSP_CLUTTER_SHADOW Two-way power factor of terrain shadowing (1 = visible).
%
%   F = rsp_clutter_shadow(L, azDeg, rangeM)
%   F = rsp_clutter_shadow(L, azDeg, rangeM, except)   ignore layer 'except'
%
%   L      : clutter layers (rsp_clutter_types)
%   azDeg  : azimuths [deg] (column), rangeM : ranges [m] (row)
%   F      : numel(azDeg) x numel(rangeM)
%
%   A layer with L.shadow (a mountain) hides what lies behind its ridge
%   (range beyond L.rangeM(2), inside its azimuth sector): the echo power
%   there is multiplied by 10^(L.shadowDb/10). The azimuth edges follow the
%   same soft ramps as the layer itself. Used for the clutter of the other
%   layers and for targets.

if nargin < 4
    except = 0;
end
F = ones(numel(azDeg), numel(rangeM));
for k = 1:numel(L)
    if ~L(k).shadow || k == except
        continue
    end
    [~, azW] = rsp_clutter_mask(L(k), azDeg, rangeM);
    behind   = rangeM(:).' > L(k).rangeM(2);
    F = F .* (1 - azW * double(behind) * (1 - 10^(L(k).shadowDb/10)));
end
end
