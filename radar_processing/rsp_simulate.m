function [video, state, truth] = rsp_simulate(P, S, pulses, state)
%RSP_SIMULATE Baseband video of a rotating search radar (block by block).
%
%   [video, state, truth] = rsp_simulate(P, S, pulses, state)
%
%   P      : chain parameters (fs, pulses, P.radar, P.antenna)
%   S      : scenario. Every level is a power PER SAMPLE at the receiver,
%            in dB on one common scale (e.g. dBm or dB re 1 unit):
%     .nRange      range samples per pulse                          (5469)
%     .noiseDb     receiver noise power                             (0 dB)
%     .blankTx     true = receiver off while transmitting
%     .seed        random seed
%     .targets     struct array, one per target:
%                    rangeM       range at time 0 [m]
%                    azDeg        azimuth [deg], 0 = north, clockwise
%                    velocityMps  radial velocity [m/s], positive = approaching
%                    powerDb      echo power at the beam peak, before pulse
%                                 compression (powerDb - noiseDb = SNR per sample)
%     .clutter     [] or struct array, one clutter layer per region:
%                    type         'ground' | 'mountain' | 'forest'
%                    powerDb, azDeg [a1 a2], rangeM [r1 r2], sigmaVMps,
%                    dcFraction, textureDb, slopeDb, shadow ... (see
%                    rsp_clutter_types for every field and the defaults)
%                  or the earlier single struct (powerDb, maxRangeM,
%                  sigmaVMps, textureDb): ground clutter at all azimuths
%            (older names still work: noisePower [linear], snrDb, cnrDb)
%   pulses : global pulse numbers of this block (1 = first pulse of the scan);
%            consecutive blocks continue the same clutter process
%   state  : [] on the first call, then the returned state
%
%   video  : numel(pulses) x nRange single complex
%   truth  : per target: azimuth, range when the beam crosses it, range
%            cell, true and folded (measurable) velocity, Doppler bin
%
%   Physics: the echo of pulse m from a target at range R(t) = R0 - v t is
%   the transmitted signal (both pulses, each at its own delay and
%   frequency) delayed by tau = 2R/c - a fractional delay applied in the
%   frequency domain - times exp(-j 2 pi f0 tau): Doppler, its sign and the
%   range walk all follow from the geometry. The amplitude follows the
%   two-way antenna pattern as the antenna turns. Noise and clutter random
%   numbers are drawn pulse by pulse, so the video does not depend on how
%   the scan is split into blocks. Intra-pulse Doppler is neglected (range-
%   Doppler coupling below 0.3 sample here). Beyond the unambiguous range
%   c/2PRF the echo of an earlier pulse arrives (second-time-around echo,
%   at delay mod PRI). Targets in the shadow of a mountain layer lose
%   shadowDb.
%
%   Clutter: each layer has a mean power map (texture x soft edges x slope
%   x shadow of mountains in front of it) smoothed by the two-way beam in
%   azimuth, and two independent parts: a static part (share dcFraction,
%   spread only by the antenna scanning) and a moving part (internal motion
%   sigmaVMps, e.g. wind in a forest). Each part is complex Gaussian noise
%   through a 4-pole slow-time filter whose bell-shaped spectrum has the
%   requested standard deviation.

c  = 299792458;
S  = rsp_levels(S);
R  = S.nRange;
N0 = 10^(S.noiseDb/10);                                         % noise power, linear
pulses = pulses(:);
M  = numel(pulses);
g  = rsp_geometry(P, R, pulses);

% --- transmitted signal: every pulse at its own delay ----------------------------
nP = numel(P.pulse);
tx = zeros(0, 1);
for k = 1:nP
    s  = rsp_waveform(P.pulse(k), P.fs);
    d  = round(P.pulse(k).delayUs * P.fs);
    n  = d + numel(s);
    tx(end+1:n, 1) = 0;
    tx(d+1:n) = tx(d+1:n) + s;
end
txLen = numel(tx);
nfft  = 2^nextpow2(R + txLen);
f     = (mod((0:nfft-1) + nfft/2, nfft) - nfft/2) * P.fs / nfft;   % MHz
TX    = fft(tx, nfft).';

% --- state (random generator, clutter process) -------------------------------------
if nargin < 4 || isempty(state)
    if isfield(S, 'seed') && ~isempty(S.seed)
        rng(S.seed);
    end
    state = clutterInit(P, S, g, R);
end
L = state.L;

SPEC = complex(zeros(M, nfft));

% random numbers drawn pulse by pulse (one column each), so the video does
% not depend on how the pulses are split into blocks
nZ = 2 * (state.nRand + R);
Z  = randn(nZ, M);

% --- clutter: every layer = static part + moving part, continuous over blocks --------
if state.cells > 0
    field = complex(zeros(M, state.cells));
    a  = g.azDeg / state.azStep;                                % texture, linear in azimuth
    i0 = floor(a);
    fr = a - i0;
    i0 = mod(i0, state.nAz);
    o  = 0;
    for k = 1:numel(state.layer)
        ly = state.layer(k);
        n  = ly.nc;
        z = zeros(M, n);
        for j = 1:2                                             % 1 static, 2 moving
            if ly.share(j) == 0
                o = o + n;
                continue
            end
            w = complex(Z(o + (1:n), :), Z(nZ/2 + o + (1:n), :)).' / sqrt(2);
            o = o + n;
            [y, state.layer(k).zi{j}] = filter(ly.b{j}, ly.a{j}, w, state.layer(k).zi{j}, 1);
            z = z + sqrt(ly.share(j)) * y;
        end
        amp = ly.amp(i0 + 1, :) .* (1 - fr) + ly.amp(i0 + 2, :) .* fr;
        cells = ly.c0 + (0:n-1);
        field(:, cells) = field(:, cells) + z .* amp;
    end
    SPEC = SPEC + fft(field / sqrt(sum(abs(tx).^2)), nfft, 2) .* TX;   % per-sample power = powerDb
end

% --- targets ----------------------------------------------------------------------------
lossLong = sum(abs(tx).^2) * P.fft.nPulses;                     % coherent gain, upper bound
priUs    = 1e6 / P.radar.prfHz;
for k = 1:numel(S.targets)
    tg  = S.targets(k);
    G   = rsp_antenna_pattern(g.azDeg - tg.azDeg, P.antenna);   % one-way power
    a   = sqrt(10^(tg.powerDb/10)) * G;                        % two-way amplitude
    if max(a)^2 * lossLong < 0.01 * N0
        continue                                                % far below the noise
    end
    Rm  = tg.rangeM - tg.velocityMps * g.timeS;                 % range at each pulse
    if ~isempty(L)
        a = a .* sqrt(rsp_clutter_shadow(L, tg.azDeg, Rm).');   % behind a mountain
    end
    tau = 2 * (Rm - P.radar.rangeOffsetM) / c * 1e6;            % delay [us]
    tauR = mod(tau, priUs);                                     % echo of an earlier pulse
                                                                % (beyond c/2PRF) lands here
    if any(tauR * P.fs + txLen > nfft - 1) || any(Rm < 0)
        continue                                                % outside the record
    end
    ph  = exp(-1j*2*pi*mod(P.radar.fcMHz * tau, 1));            % carrier phase (Doppler)
    SPEC = SPEC + (a .* ph) .* exp(-1j*2*pi*tauR * f) .* TX;
end

% --- receiver ------------------------------------------------------------------------------
rx    = ifft(SPEC, [], 2);
noise = complex(Z(state.nRand + (1:R), :), Z(nZ/2 + state.nRand + (1:R), :)).';
video = rx(:, 1:R) + sqrt(N0/2) * noise;
if S.blankTx
    video(:, 1:min(txLen, R)) = 0;
end
video = single(video);

if nargout >= 3
    truth = rsp_truth(P, S);
end
end

% ========================================================================================
function st = clutterInit(P, S, g, R)
% Clutter layers: mean power map on an azimuth x range grid (texture, soft
% edges, slope, shadow, smoothed by the two-way beam) and the slow-time
% filters of the static and moving parts.
st = struct('L', [], 'layer', [], 'cells', 0, 'nRand', 0, 'azStep', 1, 'nAz', 360);
if ~isfield(S, 'clutter') || isempty(S.clutter)
    st.L = rsp_clutter_types([]);
    return
end
L    = rsp_clutter_types(S.clutter);
st.L = L;
A    = P.antenna;
st.azStep = min(1, A.beamwidthDeg / 4);
st.nAz    = round(360 / st.azStep);
st.azStep = 360 / st.nAz;
azG  = (0:st.nAz-1)' * st.azStep;
% two-way beam (power) over the grid, for the azimuth smoothing
kern = rsp_antenna_pattern(azG, A).^2;
KERN = fft(kern / sum(kern));
% spectral spread of the antenna scanning [cycles / pulse]
sfScan = 0.265 * g.degPerPulse / A.beamwidthDeg;
vToF   = 2 / g.lambdaM / P.radar.prfHz;                         % m/s -> cycles / pulse

layer = struct('c0', {}, 'nc', {}, 'amp', {}, 'share', {}, 'b', {}, 'a', {}, 'zi', {});
for k = 1:numel(L)
    c0 = max(1, ceil((L(k).rangeM(1) - P.radar.rangeOffsetM) / g.cellM) + 1);
    c1 = min(R, floor((L(k).rangeM(2) - P.radar.rangeOffsetM) / g.cellM) + 1);
    if c1 < c0
        continue
    end
    nc = c1 - c0 + 1;
    rg = g.rangeM(c0:c1)';
    % log-normal texture, smoothed over 3 cells in range, unit mean power
    t  = randn(st.nAz, nc) * L(k).textureDb / (20*log10(exp(1)));
    t  = filter(ones(1, 3)/3, 1, t, [], 2);
    tx = exp(2*t);
    tx = tx / mean(tx(:));
    pw = 10^(L(k).powerDb/10) * rsp_clutter_mask(L(k), azG, rg) .* tx .* ...
         rsp_clutter_shadow(L, azG, rg, k);
    pw = max(real(ifft(fft(pw) .* KERN)), 0);                  % seen through the beam
    amp = sqrt(pw);
    amp(end+1, :) = amp(1, :);                                  %#ok<AGROW> wrap 360 -> 0

    % static part: only the scanning spread; moving part: internal motion too
    sfs = [0, L(k).sigmaVMps * vToF];
    if L(k).scanSpread
        sfs = sqrt(sfs.^2 + sfScan^2);
    end
    share = [L(k).dcFraction, 1 - L(k).dcFraction];
    i = numel(layer) + 1;
    layer(i).c0  = c0;
    layer(i).nc  = nc;
    layer(i).amp = amp;
    layer(i).share = share;
    layer(i).b  = cell(1, 2);
    layer(i).a  = cell(1, 2);
    layer(i).zi = cell(1, 2);
    for j = 1:2
        [layer(i).b{j}, layer(i).a{j}] = spectrumFilter(sfs(j));
        layer(i).zi{j} = [];
        if share(j) > 0                                         % start stationary
            layer(i).zi{j} = warmUp(layer(i).b{j}, layer(i).a{j}, nc);
        end
    end
    st.cells = max(st.cells, c1);
end
st.layer = layer;
st.nRand = 2 * sum([layer.nc]);
end

function [b, a] = spectrumFilter(sf)
% Four cascaded one-pole low-pass filters: a bell-shaped clutter spectrum
% with standard deviation sf [cycles / pulse] and unit power gain (the
% spectrum (1 + (f/fc)^2)^-4 has std fc / sqrt(5)).
sf = max(sf, 1e-4);
p  = exp(-2*pi * sqrt(5) * sf);
a  = 1;
for i = 1:4
    a = conv(a, [1 -p]);
end
h  = filter(1, a, [1 zeros(1, ceil(60 / (1 - p)))]);
b  = 1 / norm(h);
end

function zi = warmUp(b, a, n)
% Filter state after enough white noise to be stationary.
p  = -a(2) / 4;
nW = min(ceil(12 / (1 - p)), 20000);
zi = complex(zeros(numel(a) - 1, n));
for s = 1:500:nW
    m = min(500, nW - s + 1);
    w = complex(randn(m, n), randn(m, n)) / sqrt(2);
    [~, zi] = filter(b, a, w, zi, 1);
end
end
