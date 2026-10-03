function [value, bin] = rsp_cfar_max(map, bins)
%RSP_CFAR_MAX Output stage: maximum over the selected Doppler bins after CFAR.
%
%   [value, bin] = rsp_cfar_max(map, bins)
%
%   map   : frames x range x nfft, CFAR output (value on detections, 0 elsewhere)
%   bins  : Doppler bins taken into the output ([] = all). One bin = that
%           bin's output only, e.g. bins = 5 or bins = [2:8 10:16]
%   value : frames x range, largest detection over the selected bins (0 = none)
%   bin   : frames x range uint8, Doppler bin of that value (0 = none)

nB = size(map, 3);
if nargin < 2 || isempty(bins)
    bins = 1:nB;
end
bins = bins(:)';
if any(bins < 1 | bins > nB)
    error('rsp_cfar_max:bins', 'Bins must be within 1..%d.', nB);
end
[value, k] = max(map(:, :, bins), [], 3);
bin = uint8(bins(k));
bin(value <= 0) = 0;
end
