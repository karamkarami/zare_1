function clim = mc_clim(v)
%MC_CLIM Colour range [dB] of an image in dB (-Inf = zero).
%
%   Dense image : median .. 99.9 percentile (noise floor to the strongest echoes).
%   CFAR map (mostly zeros): every detected value, minimum .. maximum.

nz = v(v > -Inf);
if numel(nz) < 0.5 * numel(v)
    clim = [min(nz) max(nz)];
else
    clim = mc_prctile(nz, [50 99.9]);
end
if isempty(clim) || clim(2) <= clim(1)
    if isempty(clim)
        clim = 0;
    end
    clim = clim(1) + [-1 1];
end
clim = double(clim(:).');
end
