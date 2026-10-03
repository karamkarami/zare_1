function m = rsp_cfar_max(cfar, bins, S, law)
%RSP_CFAR_MAX Output stage: maximum over the selected Doppler bins after the CFAR.
%
%   m = rsp_cfar_max(cfar, bins, S, law)
%
%   cfar  : rsp_cfar output (map, threshold), or just the map (frames x range x nfft)
%   bins  : Doppler bins taken into the output ([] = all). One bin = that
%           bin's output only, e.g. bins = 5 or bins = [2:8 10:16]
%   S     : integrated values (rsp_nci output), optional
%   law   : detector law of S ('square' | 'linear' | 'log'), default 'square'
%
%   m.value     frames x range, largest CFAR output over the selected bins (0 = none)
%   m.bin       frames x range uint8, Doppler bin of that value (0 = none)
%   with S:
%   m.binFrac   frames x range single, Doppler bin refined between bins by a
%               3-point Gaussian fit on S around m.bin (NaN = none); gives the
%               velocity with a fraction of a bin accuracy
%   m.marginDb  frames x range single, S over the CFAR threshold in dB:
%               how far the detection clears the threshold (NaN = none)
%
%   This is the only maximum taken in the chain.

if isstruct(cfar)
    map = cfar.map;
else
    map = cfar;
end
nB = size(map, 3);
if nargin < 2 || isempty(bins)
    bins = 1:nB;
end
if nargin < 4 || isempty(law)
    law = 'square';
end
bins = bins(:)';
if any(bins < 1 | bins > nB)
    error('rsp_cfar_max:bins', 'Bins must be within 1..%d.', nB);
end

[m.value, k] = max(map(:, :, bins), [], 3);
on    = m.value > 0;
m.bin = uint8(bins(k));
m.bin(~on) = 0;
if nargin < 3 || isempty(S)
    return
end

% --- refined Doppler and threshold margin, on detections only ----------------------
[nF, nR, ~] = size(map);
idx = find(on);
[f, r] = ind2sub([nF nR], idx);
b   = double(m.bin(idx));
bl  = mod(b - 2, nB) + 1;                                   % Doppler is circular
br  = mod(b, nB) + 1;
at  = @(A, bb) double(A(sub2ind(size(A), f, r, bb)));
sc  = at(S, b);
if strcmpi(law, 'log')
    y = @(v) v;                                             % already logarithmic
else
    y = @(v) log(max(v, realmin));
end
yl  = y(at(S, bl));
yc  = y(sc);
yr  = y(at(S, br));
den = yl - 2*yc + yr;
d   = 0.5 * (yl - yr) ./ den;
d(den >= 0) = 0.5 * sign(yr(den >= 0) - yl(den >= 0));      % not a peak: lean to the larger side
d   = min(max(d, -0.5), 0.5);

m.binFrac = nan(nF, nR, 'single');
m.binFrac(idx) = b + d;
m.marginDb = nan(nF, nR, 'single');
if isstruct(cfar)
    th = at(cfar.threshold, b);
    if strcmpi(law, 'log')
        m.marginDb(idx) = sc - th;
    else
        m.marginDb(idx) = 10*log10(sc ./ th) * (1 + strcmpi(law, 'linear'));
    end
end
end
