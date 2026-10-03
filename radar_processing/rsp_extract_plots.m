function plots = rsp_extract_plots(det, P, nRange)
%RSP_EXTRACT_PLOTS Group detections into plots (range, azimuth, velocity).
%
%   plots = rsp_extract_plots(det, P, nRange)
%
%   det    : detections after the max stage, struct of column vectors
%            .pulse  global pulse number      .cell   range cell
%            .bin    Doppler bin (1-based)    .value  output value (linear)
%            .azDeg  (optional) antenna azimuth of the output row, e.g. from
%                    the log; otherwise computed from P.radar
%            .binFrac  (optional) refined Doppler bin (rsp_cfar_max) -> velocity
%                    between bin centres
%            .marginDb (optional) level over the CFAR threshold
%   P      : parameters (P.plots, P.radar, P.fs, P.fft)
%   nRange : range cells per pulse
%
%   plots  : struct of column vectors, one row per plot
%            .azDeg        amplitude-weighted azimuth centroid [deg]
%            .rangeM       amplitude-weighted range centroid [m]
%            .cell         range cell of the centroid
%            .bin          Doppler bin of the strongest detection
%            .velocityMps  radial velocity [m/s], positive = approaching,
%                          folded into +-lambda*PRF/4; from the refined bins
%                          (det.binFrac) when given, else the bin centre
%            .marginDb     largest level over the CFAR threshold [dB]
%            .powerDb      strongest value [dB]
%            .hits         number of detections
%            .pulse        centroid pulse number (output row, chain delay included)
%            .widthDeg     azimuth extent [deg]
%
%   Detections closer than P.plots.gapPulses pulses and P.plots.gapCells
%   cells (in any Doppler bin) are connected. Groups with fewer than
%   P.plots.minHits detections, or narrower than P.plots.minWidthDeg in
%   azimuth, are dropped: a noise alarm lasts a few pulses, a target about
%   one beamwidth. No toolbox is needed.

g  = rsp_geometry(P, nRange);
N  = numel(det.pulse);
nfft = max(P.fft.nfft, P.fft.nPulses);
empty = struct('azDeg', zeros(0, 1), 'rangeM', zeros(0, 1), 'cell', zeros(0, 1), ...
               'bin', zeros(0, 1), 'velocityMps', zeros(0, 1), 'powerDb', zeros(0, 1), ...
               'marginDb', zeros(0, 1), 'hits', zeros(0, 1), 'pulse', zeros(0, 1), ...
               'widthDeg', zeros(0, 1));
if N == 0
    plots = empty;
    return
end
p = double(det.pulse(:));
c = double(det.cell(:));

% --- neighbour pairs via a sparse lookup table --------------------------------------
p0  = min(p) - 1;
T   = sparse(p - p0, c, 1:N, max(p) - p0 + P.plots.gapPulses, nRange + P.plots.gapCells);
I = [];  J = [];
for dp = 0:P.plots.gapPulses
    for dc = -P.plots.gapCells:P.plots.gapCells
        if dp == 0 && dc <= 0
            continue                                     % each pair once
        end
        cc = c + dc;
        ok = cc >= 1;
        j  = zeros(N, 1);
        j(ok) = full(T(sub2ind(size(T), p(ok) - p0 + dp, cc(ok))));
        k  = find(j > 0);
        I  = [I; k];                                     %#ok<AGROW>
        J  = [J; j(k)];                                  %#ok<AGROW>
    end
end

% --- connected components: minimum-label propagation with pointer jumping -----------
lab = (1:N)';
while true
    old = lab;
    m   = accumarray([I; J], [lab(J); lab(I)], [N 1], @min, Inf);
    lab = min(lab, m);
    lab = lab(lab);
    lab = lab(lab);
    if isequal(lab, old)
        break
    end
end
[~, ~, grp] = unique(lab);

% --- plot parameters -----------------------------------------------------------------
hits = accumarray(grp, 1);
span = (accumarray(grp, p, [], @max) - accumarray(grp, p, [], @min) + 1) * g.degPerPulse;
keep = find(hits >= P.plots.minHits & span >= P.plots.minWidthDeg);
nPl  = numel(keep);
plots = empty;
w   = double(det.value(:));
% azimuth of the data centre behind each output row (chain delay removed)
if isfield(det, 'azDeg') && ~isempty(det.azDeg)
    az = mod(double(det.azDeg(:)) - g.degPerPulse * g.outputDelayPulses, 360);
else
    az = mod(P.radar.azStartDeg + g.degPerPulse * (p - 1 - g.outputDelayPulses), 360);
end
for i = 1:nPl
    k  = find(grp == keep(i));
    wk = w(k);
    a0 = az(k(1));
    da = mod(az(k) - a0 + 180, 360) - 180;              % unwrap around the first hit
    [pk, im] = max(wk);
    plots.azDeg(i, 1)       = mod(a0 + sum(wk .* da) / sum(wk), 360);
    plots.cell(i, 1)        = sum(wk .* c(k)) / sum(wk);
    plots.rangeM(i, 1)      = P.radar.rangeOffsetM + (plots.cell(i) - 1) * g.cellM;
    plots.bin(i, 1)         = double(det.bin(k(im)));
    if isfield(det, 'binFrac') && ~isempty(det.binFrac)
        % power-weighted circular mean of the refined bins (detections within
        % 10 dB of the strongest one)
        sel = wk >= pk / 10;
        ph  = 2*pi * (double(det.binFrac(k(sel))) - 1) / nfft;
        fb  = angle(sum(wk(sel) .* exp(1j*ph))) / (2*pi);      % cycles / pulse
        fb  = fb - P.fft.shift * floor(nfft/2) / nfft;
        fb  = mod(fb + 0.5, 1) - 0.5;
        plots.velocityMps(i, 1) = fb * P.radar.prfHz * g.lambdaM / 2;
    else
        plots.velocityMps(i, 1) = g.binVelMps(plots.bin(i));
    end
    if isfield(det, 'marginDb') && ~isempty(det.marginDb)
        plots.marginDb(i, 1) = max(double(det.marginDb(k)));
    else
        plots.marginDb(i, 1) = NaN;
    end
    plots.powerDb(i, 1)     = 10*log10(pk);
    plots.hits(i, 1)        = numel(k);
    plots.pulse(i, 1)       = sum(wk .* p(k)) / sum(wk);
    plots.widthDeg(i, 1)    = max(da) - min(da) + g.degPerPulse;
end
end
