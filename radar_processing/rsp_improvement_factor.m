function IF = rsp_improvement_factor(data, P, varargin)
%RSP_IMPROVEMENT_FACTOR MTI improvement factor measured on recorded clutter.
%
%   IF = rsp_improvement_factor(video, P)
%   IF = rsp_improvement_factor(video, P, 'option', value, ...)
%
%   video : consecutive pulses x range samples, complex baseband. Raw video of
%           a few hundred PRIs is enough, and the antenna may turn.
%   P     : parameters (radar_params): the two pulses, PRF, canceler, FFT
%
%   Definition (IEEE): the signal-to-clutter ratio at the output of the
%   clutter filter over the one at its input, averaged over all target
%   velocities. Averaged over velocity, a target gets the filter's noise
%   gain G = sum|w|^2, so
%
%       I = G * Cin / Cout = CNRin / CNRout
%
%   Cin = Pin - N is the clutter power at the input and Cout = Pout - N*G the
%   clutter residue at the output. Only the noise power N is needed, and the
%   clutter does not have to be located by hand:
%     1. noise cells   : no pulse-to-pulse correlation (|rho(1)| < 3/sqrt(pulses)).
%                        N is their median power.
%     2. clutter cells : CNR >= minCnrDb and mean radial velocity inside
%                        +-maxVelocityMps (stationary clutter). Cells where a
%                        moving echo passes the canceler or a Doppler bin are
%                        left out, with their neighbours within the range
%                        resolution: a block of 32 pulses whose residue is
%                        10 dB above what the cell's clutter gives.
%     3. I             : pooled over the clutter cells (summed powers), per
%                        pulse zone (short / long pulse) and per range cell.
%                        Also for 2/3/4-pulse cancelers and for every Doppler bin.
%   Measurement floor: the residue is measured against the noise, and the
%   noise estimate has a spread that falls with the number of pulses and cells.
%   A residue below 3 standard deviations cannot be told from zero. I is then
%   only a lower bound (bound = true), and maxDb is the highest I this data
%   can show.
%
%   Options
%     'input'          'video'  'video'   : raw video, decoded here (both pulses,
%                                           rsp_matched_filter, stitched)
%                               'decoded' : already the decoder output (log lane)
%                               'samples' : raw samples, not decoded (codes not
%                                           known; less CNR, so a lower maxDb)
%     'pulses'         []       rows used ([] = all)
%     'cells'          []       [first last] range cells ([] = all valid cells)
%     'noise'          []       noise power per cell at the canceler input, one
%                               value or one per range cell ([] = estimated)
%     'noiseCells'     []       [first last] cells holding noise only ([] = found)
%     'minCnrDb'       10       clutter cells: CNR at least this
%     'maxVelocityMps' 5        clutter cells: |mean radial velocity| at most this
%     'orders'         [2 3 4]  canceler orders compared on the same clutter
%     'plot'           true
%     'quiet'          false    true = no printing
%
%   IF fields
%     .cells, .rangeM   analysed range cells and their range [m]
%     .noise            noise power per cell; .noiseCells (logical)
%     .noiseInfo        .method .nCells .rangeCorrelation (noise correlation
%                       factor along range) .relStd (of the noise estimate)
%     .cell             per cell: .pIn .pOut (canceler output) .cnrDb .rho1
%                       .velMps .clutter (logical) .iDb .bound
%     .clutter          .nCells .rangeKm .cnrDb (pooled) .rho1 .spreadMps
%                       (Gaussian spectrum std) .scanSpreadMps (antenna scan
%                       alone) .hits (pulses per 3 dB beamwidth) .meanVelMps
%                       .nMoving (cells left out: moving echo)
%     .canceler         P.canceler: .taps .gainDb .iDb .bound .maxDb .errDb
%                       .caDb (clutter attenuation Cin/Cout) .scanLimitDb;
%                       per cell: .medianDb .p10Db .p90Db (cells that measure
%                       their I) .boundShare (cells that give only a bound);
%                       .zone(k): .name .cells .iDb .bound .maxDb .errDb .nCells
%     .orders(k)        .order .taps .iDb .bound .maxDb .errDb .scanLimitDb
%     .bins             canceler + Doppler FFT, per bin: .velMps .gainDb (average
%                       over velocity) .peakGainDb .iDb (average over velocity)
%                       .iPeakDb (target at the bin's peak response: includes
%                       the coherent gain) .bound .maxDb .taps (taps x bins)
%     .spectrum         mean Doppler spectrum [dB re noise]: .velMps .clutterDb
%                       (clutter cells) .residueDb (their clutter after the
%                       canceler, / G) .noiseDb (noise cells)
%     .mf               decoder info (input 'video')

opt = struct('input', 'video', 'pulses', [], 'cells', [], 'noise', [], 'noiseCells', [], ...
             'minCnrDb', 10, 'maxVelocityMps', 5, 'orders', [2 3 4], 'plot', true, ...
             'quiet', false);
for i = 1:2:numel(varargin)
    if ~ischar(varargin{i}) || ~isfield(opt, varargin{i})
        error('rsp_improvement_factor:option', 'Unknown option number %d.', (i + 1)/2);
    end
    opt.(varargin{i}) = varargin{i+1};
end

if ~isempty(opt.pulses)
    data = data(opt.pulses, :);
end
if ~isfloat(data)
    data = single(data);
end
[M, R] = size(data);
geo   = rsp_geometry(P, R);
toVel = geo.lambdaM / 2 * P.radar.prfHz;                % m/s per cycle/pulse

C = P.canceler;
C.output = 'valid';
w = rsp_canceler_taps(C);
if ~C.enable
    w = 1;
end
G    = sum(abs(w).^2);                                  % noise gain = mean signal gain
nMin = numel(w) + P.fft.nPulses + 16;
if M < nMin
    error('rsp_improvement_factor:pulses', 'Need at least %d consecutive pulses, got %d.', nMin, M);
end

%% Decoder ---------------------------------------------------------------------------
mf = [];
switch lower(opt.input)
    case 'video'
        [X, mf] = rsp_matched_filter(data, P);
        valid   = mf.validCells;
    case {'decoded', 'samples'}
        X     = data;
        valid = [1 R];
    otherwise
        error('rsp_improvement_factor:input', 'Unknown input "%s".', opt.input);
end
clear data
if ~isempty(opt.cells)
    valid = [max(valid(1), opt.cells(1)) min(valid(2), opt.cells(end))];
end
cells = valid(1):valid(2);
X     = double(X(:, cells));
nC    = numel(cells);

%% Per cell: power, pulse-to-pulse correlation, mean Doppler ----------------------------
pIn  = mean(abs(X).^2, 1);
r1   = mean(X(2:end, :) .* conj(X(1:end-1, :)), 1);
rho1 = abs(r1) ./ max(pIn, realmin);
vel  = angle(r1) / (2*pi) * toVel;                      % positive = approaching
live = pIn > 0;                                         % blanked cells are all zero

%% Noise -------------------------------------------------------------------------------
if isempty(opt.noiseCells)
    cand = live & rho1 < 3 / sqrt(M - 1);               % white along pulses
    how  = 'cells without pulse-to-pulse correlation';
else
    cand = live & cells >= opt.noiseCells(1) & cells <= opt.noiseCells(end);
    how  = sprintf('cells %d-%d', opt.noiseCells(1), opt.noiseCells(end));
end
if isempty(opt.noise)
    N = repmat(prc(pIn(cand), 50), 1, nC);
else
    N = double(opt.noise(:)).';
    if isscalar(N)
        N = repmat(N, 1, nC);
    elseif numel(N) == R
        N = N(cells);
    else
        error('rsp_improvement_factor:noise', 'Give one noise power or one per range cell (%d).', R);
    end
    how = 'given';
end
isNoise = cand & abs(pIn - N) < 6 * N / sqrt(M);       % drops weak clutter and echoes
nNoise  = nnz(isNoise);
if isempty(opt.noise)
    if nNoise < 10
        error('rsp_improvement_factor:noise', ['Only %d noise-only cells found. Give the ' ...
              'noise power (''noise'') or a clutter-free range (''noiseCells'').'], nNoise);
    end
    N(:) = prc(pIn(isNoise), 50);
end
kappa = rangeCorrelation(X(:, isNoise) ./ sqrt(N(isNoise)), find(isNoise));
if isempty(opt.noise)
    relVarN = pi/2 * kappa / (M * nNoise);              % variance of a median
else
    relVarN = 0;
end

%% Canceler, Doppler bins and clutter cells -------------------------------------------------
Y     = rsp_canceler(X, C);                             % rows: pulses numel(w)..M
F     = size(Y, 1);
pOut  = mean(abs(Y).^2, 1);
phi   = noiseCorrelation(w, 1);
cnr   = (pIn - N) ./ N;
S     = live & ~isNoise & cnr >= 10^(opt.minCnrDb/10); % strong echo ...
still = abs(vel) <= opt.maxVelocityMps;                % ... standing still
Xa    = abs(X(end-F+1:end, :)).^2;                     % input power, rows of Y
moving = movingEcho(Xa, abs(Y).^2, N, N * G, phi);

% every Doppler bin of these cells: a moving echo near the canceler notch is
% weak after the canceler but strong in its own bin
[h, gainB, peakB] = binFilters(P, C);
nB   = size(h, 2);
phiB = zeros(1, nB);
for j = 1:nB
    phiB(j) = noiseCorrelation(h(:, j), P.fft.hop);
end
iS = find(S & still);
pB = zeros(numel(iS), nB);
FB = 0;
for c0 = 1:256:numel(iS)
    j = c0:min(c0 + 255, numel(iS));
    q = iS(j);
    [D, iD] = rsp_doppler_fft(Y(:, q), P.fft);
    FB = size(D, 1);
    pD = abs(D).^2;
    pB(j, :) = reshape(mean(pD, 1), numel(j), nB);
    for k = 1:nB
        moving(q) = moving(q) | movingEcho(Xa(iD, q), pD(:, :, k), N(q), ...
                                           N(q) * gainB(k), phiB(k));
    end
end
clear D pD
K       = ceil(kappa);                                  % the echo's neighbours in range
moving  = conv(double(moving), ones(1, 2*K + 1), 'same') > 0;
nMoving = nnz(S & ~(still & ~moving));
S       = S & still & ~moving;
nS      = nnz(S);
pB      = pB(S(iS), :);                                 % bins of the clutter cells

nHits = P.antenna.beamwidthDeg / geo.degPerPulse;      % pulses per 3 dB beamwidth
if geo.degPerPulse == 0
    nHits = Inf;
end

IF.cells  = cells;
IF.rangeM = geo.rangeM(cells).';
IF.noise  = N;
IF.noiseCells = isNoise;
IF.noiseInfo  = struct('method', how, 'nCells', nNoise, 'rangeCorrelation', kappa, ...
                       'relStd', sqrt(relVarN));
IF.cell = struct('pIn', pIn, 'pOut', pOut, 'cnrDb', 10*log10(max(cnr, realmin)), ...
                 'rho1', rho1, 'velMps', vel, 'clutter', S, 'iDb', nan(1, nC), ...
                 'bound', false(1, nC));
IF.mf = mf;

cin = pIn - N;
c.nCells        = nS;
c.rangeKm       = [NaN NaN];
c.cnrDb         = NaN;
c.rho1          = NaN;
c.spreadMps     = NaN;
c.scanSpreadMps = toVel * sqrt(log(2)) / (pi * nHits);
c.hits          = nHits;
c.meanVelMps    = NaN;
c.nMoving       = nMoving;
if nS > 0
    c.rangeKm    = IF.rangeM([find(S, 1) find(S, 1, 'last')]) / 1e3;
    c.cnrDb      = 10*log10(sum(cin(S)) / sum(N(S)));
    c.rho1       = abs(sum(r1(S))) / sum(cin(S));      % clutter only (noise removed)
    c.spreadMps  = toVel * sqrt(max(-log(min(c.rho1, 1)), 0) / (2*pi^2));
    c.meanVelMps = angle(sum(r1(S))) / (2*pi) * toVel;
end
IF.clutter = c;

%% Improvement factor of the canceler ------------------------------------------------------
res     = pOut - N * G;
sigCell = N * G .* sqrt(phi / F + relVarN);
[iDb, bnd] = measure(cin, res, G, sigCell);
IF.cell.iDb(S)   = iDb(S);
IF.cell.bound(S) = bnd(S);

k = pool(S, cin, res, N * G, G, phi, F, kappa, relVarN);
k.taps        = w;
k.gainDb      = 10*log10(G);
k.scanLimitDb = scanLimit(w, nHits);
v = iDb(S & ~bnd);                                      % cells that measure their I
k.medianDb   = prc(v, 50);
k.p10Db      = prc(v, 10);
k.p90Db      = prc(v, 90);
k.boundShare = mean(bnd(S));

% zones: the short pulse covers the blind zone of the long one
zones = struct('name', {}, 'mask', {});
if ~isempty(mf) && P.mf.enable && numel(P.pulse) > 1 && strcmpi(P.mf.combine, 'stitch')
    zones(1).name = sprintf('pulse %d (%s)', mf.shortPulse, P.pulse(mf.shortPulse).type);
    zones(1).mask = cells <= mf.switchCell;
    zones(2).name = sprintf('pulse %d (%s)', mf.longPulse, P.pulse(mf.longPulse).type);
    zones(2).mask = cells > mf.switchCell;
else
    zones(1).name = 'all cells';
    zones(1).mask = true(1, nC);
end
k.zone = struct('name', {}, 'cells', {}, 'iDb', {}, 'bound', {}, 'maxDb', {}, ...
                'errDb', {}, 'nCells', {});
zones = zones(arrayfun(@(z) any(z.mask), zones));      % zones inside the analysed cells
for z = 1:numel(zones)
    r = pool(S & zones(z).mask, cin, res, N * G, G, phi, F, kappa, relVarN);
    k.zone(end+1) = struct('name', zones(z).name, ...
                           'cells', cells([find(zones(z).mask, 1) find(zones(z).mask, 1, 'last')]), ...
                           'iDb', r.iDb, 'bound', r.bound, 'maxDb', r.maxDb, ...
                           'errDb', r.errDb, 'nCells', r.nCells);
end
IF.canceler = k;

%% Other canceler orders on the same clutter ---------------------------------------------
IF.orders = struct('order', {}, 'taps', {}, 'iDb', {}, 'bound', {}, 'maxDb', {}, ...
                   'errDb', {}, 'scanLimitDb', {});
IF.bins   = struct();
IF.spectrum = struct();
if nS == 0
    if ~opt.quiet
        fprintf(['\nImprovement factor: no clutter cells (CNR >= %g dB, |v| <= %g m/s). ' ...
                 'Lower ''minCnrDb'' or check the noise (%.1f dB).\n'], opt.minCnrDb, ...
                opt.maxVelocityMps, 10*log10(median(N)));
    end
    if opt.plot
        plotIF(IF, P, zones);
    end
    return
end
XS = X(:, S);
NS = N(S);
for j = 1:numel(opt.orders)
    Ck = C;
    Ck.enable    = true;
    Ck.order     = opt.orders(j);
    Ck.weights   = [];
    Ck.normalize = false;
    wk = rsp_canceler_taps(Ck);
    Yk = rsp_canceler(XS, Ck);
    Gk = sum(wk.^2);
    r  = pool(true(1, nS), cin(S), mean(abs(Yk).^2, 1) - NS * Gk, NS * Gk, Gk, ...
              noiseCorrelation(wk, 1), size(Yk, 1), kappa, relVarN);
    IF.orders(j) = struct('order', opt.orders(j), 'taps', wk, 'iDb', r.iDb, ...
                          'bound', r.bound, 'maxDb', r.maxDb, 'errDb', r.errDb, ...
                          'scanLimitDb', scanLimit(wk, nHits));
end

%% Canceler + Doppler FFT, every bin ----------------------------------------------------------
b = struct('velMps', geo.binVelMps, 'gainDb', 10*log10(gainB), 'peakGainDb', 10*log10(peakB), ...
           'iDb', zeros(1, nB), 'iPeakDb', zeros(1, nB), 'bound', false(1, nB), ...
           'maxDb', zeros(1, nB), 'taps', h);
for j = 1:nB
    r = pool(true(1, nS), cin(S), pB(:, j).' - NS * gainB(j), NS * gainB(j), gainB(j), ...
             phiB(j), FB, kappa, relVarN);
    b.iDb(j)     = r.iDb;
    b.iPeakDb(j) = r.iDb + 10*log10(peakB(j) / gainB(j));
    b.bound(j)   = r.bound;
    b.maxDb(j)   = r.maxDb;
end
IF.bins = b;

%% Mean Doppler spectrum (dB re noise) --------------------------------------------------------
[f, sC] = meanSpectrum(XS ./ sqrt(NS));
[~, sN] = meanSpectrum(X(:, isNoise) ./ sqrt(N(isNoise)));
W = fftshift(abs(fft(w(:), numel(f))).^2) / G;
IF.spectrum = struct('velMps', f * toVel, 'clutterDb', 10*log10(sC), ...
                     'residueDb', 10*log10(max(sC - 1, 0) .* W + realmin), ...
                     'noiseDb', 10*log10(sN));

if ~opt.quiet
    report(IF, opt, P, M, mf);
end
if opt.plot
    plotIF(IF, P, zones);
end
end

% ======================================================================================
function [iDb, bound, maxDb, errDb] = measure(cin, res, G, sigma)
% I = G*cin/res, or a lower bound when the residue is below the floor (3 sigma).
fl    = 3 * sigma;
bound = res < fl;
iDb   = 10*log10(G * cin ./ max(res, fl));
maxDb = 10*log10(G * cin ./ fl);
errDb = 5*log10((res + sigma) ./ max(res - sigma, realmin));
errDb(bound) = NaN;
end

function r = pool(S, cin, res, nOut, G, phi, frames, kappa, relVarN)
% Improvement factor of the clutter cells S together (ratio of the summed powers).
% Spread of the summed residue: noise of the filter output (phi along pulses,
% kappa along range) and the error of the noise estimate (common to all cells).
r.nCells = nnz(S);
if r.nCells == 0
    [r.iDb, r.maxDb, r.errDb, r.caDb] = deal(NaN);
    r.bound = false;
    return
end
sigma = sqrt(phi / frames * kappa * sum(nOut(S).^2) + relVarN * sum(nOut(S))^2);
[r.iDb, r.bound, r.maxDb, r.errDb] = measure(sum(cin(S)), sum(res(S)), G, sigma);
r.caDb = r.iDb - 10*log10(G);
end

function phi = noiseCorrelation(h, hop)
% White noise through taps h, outputs hop pulses apart: sum_l |r(l)|^2 / r(0)^2,
% the variance factor of a mean output power.
h    = h(:);
L    = numel(h);
r    = conv(h, conj(flipud(h)));                         % lags -(L-1) .. L-1
keep = mod(-(L-1):(L-1), hop) == 0;
phi  = sum(abs(r(keep)).^2) / abs(r(L))^2;
end

function kappa = rangeCorrelation(Xn, cellIdx)
% Noise correlation factor along range, 1 + 2 sum_k |rho(k)|^2 (unit-power noise,
% e.g. from the decoder): the variance of a sum over cells grows by kappa.
kappa = 1;
small = 0;
for k = 1:64
    [~, a, b] = intersect(cellIdx, cellIdx + k);        % pairs of noise cells k apart
    if numel(a) < 10
        break
    end
    rho = mean(mean(Xn(:, a) .* conj(Xn(:, b))));
    if abs(rho) < 0.05
        small = small + 1;
        if small >= 3
            break
        end
    else
        small = 0;
    end
    kappa = kappa + 2 * abs(rho)^2;
end
end

function moving = movingEcho(pa, pb, N, nOut, phi)
% Cells where a moving echo passes a clutter filter. pa, pb: power at the input
% and at the output (rows aligned in time x cells), nOut: output noise power.
% A block of 32 rows whose residue is 10 dB above what the cell gives in the
% other blocks (its median residue, or its median residue / clutter ratio times
% the block's clutter power) plus the noise spread. A strong fixed echo keeps
% its ratio, so it is not flagged.
Lb = 32;
[F, nC] = size(pb);
B = floor(F / Lb);
moving = false(1, nC);
if B < 4 || nC == 0
    return
end
rows = 1:B*Lb;
a  = reshape(mean(reshape(pa(rows, :), Lb, B, nC), 1), B, nC) - N;       % clutter
b  = reshape(mean(reshape(pb(rows, :), Lb, B, nC), 1), B, nC) - nOut;    % residue
sb = nOut * sqrt(phi / Lb);
q  = max(median(b ./ max(a, realmin), 1), 0);
moving = any(b > 10 * (max(a .* q, median(b, 1)) + sb), 1);
end

function [h, gain, peak] = binFilters(P, C)
% Slow-time taps of canceler + Doppler FFT, from an impulse through rsp_canceler
% and rsp_doppler_fft: bin k of the FFT ending at pulse m is sum_i h(i+1,k) x(m-i).
F = P.fft;
F.hop = 1;
nT = (numel(rsp_canceler_taps(C)) - 1) * C.enable;
L  = F.nPulses + nT;
x  = zeros(2*L + 1, 1);
x(L + 1) = 1;                                           % after the canceler transient
[y, iy] = rsp_canceler(x, C, (1:numel(x))');
[D, iD] = rsp_doppler_fft(y, F, iy);
D = reshape(D, size(D, 1), []);
[~, row] = ismember(L + 1 + (0:L-1)', iD);
h    = D(row, :);
gain = sum(abs(h).^2, 1);                              % average over velocity
peak = max(abs(fft(h, 4096, 1)).^2, [], 1);           % at the bin's best velocity
end

function iDb = scanLimit(w, nHits)
% Limit set by the antenna scan alone (Gaussian beam, nHits pulses in the one-way
% 3 dB beamwidth): rho(k) = exp(-2 ln2 (k/nHits)^2), I = w w' / (w T w').
if ~isfinite(nHits)
    iDb = Inf;
    return
end
w = w(:).';
k = 0:numel(w)-1;
T = toeplitz(expm1(-2*log(2) * (k / nHits).^2));      % rho - 1 keeps the precision
q = real(w * T * w') + abs(sum(w))^2;
iDb = 10*log10(real(w * w') / max(q, realmin));
end

function [f, s] = meanSpectrum(Z)
% Mean Doppler spectrum (Welch, Blackman-Harris, white noise = 1 per bin).
[M, n] = size(Z);
Ls  = min(128, M);
nf  = max(256, Ls);
win = rsp_window('blackmanharris', Ls);
win = win / norm(win);
starts = 1:floor(Ls/2):(M - Ls + 1);
f = (-nf/2:nf/2-1)' / nf;                              % cycles per pulse
s = zeros(nf, 1);
for c0 = 1:256:n
    Zc = Z(:, c0:min(c0 + 255, n));
    for s0 = starts
        s = s + sum(abs(fft(Zc(s0:s0+Ls-1, :) .* win, nf, 1)).^2, 2);
    end
end
s = fftshift(s) / max(numel(starts) * n, 1);
end

function v = prc(x, p)
% Percentile without a toolbox (NaN when x is empty).
x = sort(x(isfinite(x)));
if isempty(x)
    v = NaN;
    return
end
v = x(min(max(round(p/100 * numel(x)), 1), numel(x)));
end

function s = fmtI(iDb, bound, errDb)
if bound
    s = sprintf('> %.1f dB', iDb);
elseif isnan(iDb)
    s = '-';
else
    s = sprintf('= %.1f dB (+-%.1f)', iDb, errDb);
end
end

function report(IF, opt, P, M, mf)
km = IF.rangeM([1 end]) / 1e3;
switch lower(opt.input)
    case 'video'
        src = sprintf('decoded here, %d pulses, switch cell %d', numel(P.pulse), mf.switchCell);
    otherwise
        src = opt.input;
end
c = IF.clutter;
k = IF.canceler;
fprintf('\nImprovement factor\n');
fprintf('  data     : %d pulses x %d cells (%s), %.2f-%.2f km\n', M, numel(IF.cells), src, km);
fprintf('  noise    : %.1f dB per cell (%s: %d cells, range correlation %.1f)\n', ...
        10*log10(median(IF.noise)), IF.noiseInfo.method, IF.noiseInfo.nCells, ...
        IF.noiseInfo.rangeCorrelation);
fprintf('  clutter  : %d cells, %.2f-%.2f km (CNR >= %g dB, |v| <= %g m/s; %d cells with a moving echo left out)\n', ...
        c.nCells, c.rangeKm, opt.minCnrDb, opt.maxVelocityMps, c.nMoving);
fprintf('             CNR %.1f dB, rho(1) %.5f, spread %.2f m/s (antenna scan alone %.2f m/s), mean %+.2f m/s\n', ...
        c.cnrDb, c.rho1, c.spreadMps, c.scanSpreadMps, c.meanVelMps);
fprintf('  canceler %s: I %s, clutter attenuation %.1f dB, noise gain %.1f dB\n', ...
        mat2str(k.taps, 4), fmtI(k.iDb, k.bound, k.errDb), k.caDb, k.gainDb);
fprintf('             this data can show I up to %.1f dB; antenna scan limit %.1f dB\n', ...
        k.maxDb, k.scanLimitDb);
for z = 1:numel(k.zone)
    q = k.zone(z);
    if q.nCells == 0
        fprintf('     %-14s cells %5d-%5d : no clutter cells\n', q.name, q.cells);
    elseif q.bound
        fprintf('     %-14s cells %5d-%5d : I %s, %d clutter cells (residue below the floor)\n', ...
                q.name, q.cells, fmtI(q.iDb, q.bound, q.errDb), q.nCells);
    else
        fprintf('     %-14s cells %5d-%5d : I %s, %d clutter cells (can show up to %.1f dB)\n', ...
                q.name, q.cells, fmtI(q.iDb, q.bound, q.errDb), q.nCells, q.maxDb);
    end
end
if isnan(k.medianDb)
    fprintf('     per cell: every cell is too weak and gives only a lower bound\n');
else
    fprintf(['     per cell: median %.1f dB, 10-90 %% %.1f-%.1f dB; %.0f %% of the cells are too weak ' ...
             'and give only a lower bound\n'], k.medianDb, k.p10Db, k.p90Db, 100 * k.boundShare);
end
fprintf('  canceler orders on the same clutter:\n     order  I                      scan limit\n');
for j = 1:numel(IF.orders)
    o = IF.orders(j);
    fprintf('     %3d    %-22s %6.1f dB\n', o.order, fmtI(o.iDb, o.bound, o.errDb), o.scanLimitDb);
end
b = IF.bins;
fprintf('  canceler + %d-pulse FFT (%s), per bin:\n', P.fft.nPulses, P.fft.window);
fprintf('     bin  velocity   I average     I at bin peak   coherent gain\n');
for j = 1:numel(b.iDb)
    mark = ' ';
    if b.bound(j)
        mark = '>';
    end
    fprintf('     %3d  %7.1f m/s  %s%6.1f dB    %s%6.1f dB      %5.1f dB\n', j, b.velMps(j), ...
            mark, b.iDb(j), mark, b.iPeakDb(j), b.peakGainDb(j) - b.gainDb(j));
end
fprintf('     (> = lower bound: the residue is below what this data can measure)\n');
end

function plotIF(IF, P, zones)
km = IF.rangeM / 1e3;
S  = IF.cell.clutter;
G  = 10^(IF.canceler.gainDb/10);
figure('Name', 'Improvement factor', 'Color', 'w', 'Position', [80 80 1300 800]);

subplot(2, 2, 1);
plot(km, 10*log10(IF.cell.pIn + realmin), 'Color', [0.6 0.6 0.6]);
hold on;
plot(km, 10*log10(IF.cell.pOut / G + realmin), 'b');
plot(km(S), 10*log10(IF.cell.pIn(S)), 'r.', 'MarkerSize', 6);
plot(km, 10*log10(IF.noise), 'k--', 'LineWidth', 1);
lim = [10*log10(min(IF.noise)) - 10, 10*log10(max(IF.cell.pIn)) + 5];
ylim(lim);
if numel(zones) > 1
    sw = km(find(zones(1).mask, 1, 'last'));
    plot([sw sw], lim, 'k:');
end
grid on;
xlabel('range [km]');
ylabel('power per cell [dB]');
legend({'input', 'after canceler / G', 'clutter cells', 'noise'}, 'Location', 'northeast');
title(sprintf('%d clutter cells, CNR %.1f dB', IF.clutter.nCells, IF.clutter.cnrDb));
if ~any(S)
    return
end

subplot(2, 2, 2);
m = S & ~IF.cell.bound;
b = S & IF.cell.bound;
plot(km(m), IF.cell.iDb(m), 'b.', 'MarkerSize', 6);
hold on;
plot(km(b), IF.cell.iDb(b), '^', 'Color', [0.3 0.75 0.9], 'MarkerSize', 3);
for z = 1:numel(zones)                                 % pooled I of each pulse zone
    q = find(S & zones(z).mask);
    if ~isempty(q)
        style = 'r-';
        if IF.canceler.zone(z).bound
            style = 'r--';                              % a lower bound
        end
        plot(km(q([1 end])), IF.canceler.zone(z).iDb * [1 1], style, 'LineWidth', 2);
    end
end
grid on;
xlabel('range [km]');
ylabel('I [dB]');
legend({'per cell', 'per cell, lower bound', 'pulse zone (dashed = bound)'}, ...
       'Location', 'southwest');
title(sprintf('Canceler %s: I %s', mat2str(IF.canceler.taps, 3), ...
      fmtI(IF.canceler.iDb, IF.canceler.bound, IF.canceler.errDb)));

subplot(2, 2, 3);
plot(IF.spectrum.velMps, IF.spectrum.clutterDb, 'r');
hold on;
plot(IF.spectrum.velMps, IF.spectrum.residueDb, 'b');
plot(IF.spectrum.velMps, IF.spectrum.noiseDb, 'Color', [0.5 0.5 0.5]);
grid on;
xlim(IF.spectrum.velMps([1 end]));
ylim([-20, max(IF.spectrum.clutterDb) + 5]);
xlabel('radial velocity [m/s]');
ylabel('dB re noise');
legend({'clutter cells', 'their clutter after canceler / G', 'noise cells'}, 'Location', 'northwest');
title(sprintf('Mean Doppler spectrum (spread %.2f m/s)', IF.clutter.spreadMps));

subplot(2, 2, 4);
[v, o] = sort(IF.bins.velMps);
plot(v, IF.bins.iDb(o), 'bo-');
hold on;
plot(v, IF.bins.iPeakDb(o), 'rs-');
plot(v([1 end]), IF.canceler.iDb * [1 1], 'k--');
bd = IF.bins.bound(o);
plot(v(bd), IF.bins.iDb(o(bd)), 'b^', 'MarkerFaceColor', 'b');
grid on;
xlabel('Doppler bin velocity [m/s]');
ylabel('I [dB]');
legend({'average over velocity', 'target at the bin peak', 'canceler alone'}, ...
       'Location', 'southwest');
title(sprintf('Canceler + %d-pulse FFT (%s), filled = lower bound', P.fft.nPulses, P.fft.window));
end
