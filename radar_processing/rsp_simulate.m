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
%     .clutter     [] or struct (ground clutter, all azimuths):
%                    powerDb      mean clutter power
%                    maxRangeM    clutter extent (power fades over the last 30 %)
%                    sigmaVMps    spectral spread (wind + scanning) [m/s]
%                    textureDb    spatial fluctuation (log-normal) [dB]
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
%   Doppler coupling below 0.3 sample here).

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
    state = struct('zi', [], 'texture', [], 'cells', 0, 'h', []);
    if isfield(S, 'clutter') && ~isempty(S.clutter)
        cl = S.clutter;
        state.cells   = min(R, floor((cl.maxRangeM - P.radar.rangeOffsetM) / g.cellM) + 1);
        % log-normal texture on a 1 deg x cell grid, smoothed over 3 cells
        t = randn(361, state.cells) * cl.textureDb / (20*log10(exp(1)));
        t = filter(ones(1, 3)/3, 1, t, [], 2);
        t(361, :) = t(1, :);                                   % wrap 360 -> 0
        state.texture = exp(t);
        state.texture = state.texture / sqrt(mean(state.texture(:).^2));   % unit mean power
        % fade out over the last 30 % of the extent (no hard clutter edge)
        r = (0:state.cells-1) / max(state.cells - 1, 1);
        state.texture = state.texture .* min(1, cos(pi/2 * max(r - 0.7, 0) / 0.3).^2);
        % Gaussian clutter spectrum (std sf cycles/pulse) = Gaussian FIR in slow time
        sf  = 2 * cl.sigmaVMps / g.lambdaM / P.radar.prfHz;
        st  = max(1 / (2*pi*sf*sqrt(2)), 0.5);                  % impulse response std [pulses]
        m   = (-ceil(4*st):ceil(4*st))';
        hc  = exp(-m.^2 / (2*st^2));
        state.h  = hc / norm(hc);                               % unit power gain
        state.zi = complex(zeros(numel(hc) - 1, state.cells));
        % warm-up: run the filter once so the first block is already stationary
        w = (randn(numel(hc), state.cells) + 1j*randn(numel(hc), state.cells)) / sqrt(2);
        [~, state.zi] = filter(state.h, 1, w, state.zi, 1);
    end
end

SPEC = complex(zeros(M, nfft));

% random numbers drawn pulse by pulse (one column each), so the video does
% not depend on how the pulses are split into blocks
nZ = 2 * (state.cells + R);
Z  = randn(nZ, M);

% --- clutter: zero-mean Doppler, Gaussian spectrum, continuous over blocks ------------
if state.cells > 0
    cl  = S.clutter;
    w   = complex(Z(1:state.cells, :), Z(nZ/2 + (1:state.cells), :)).' / sqrt(2);
    [z, state.zi] = filter(state.h, 1, w, state.zi, 1);
    a0  = floor(g.azDeg);                                       % texture, bilinear in azimuth
    fr  = g.azDeg - a0;
    tex = state.texture(a0 + 1, :) .* (1 - fr) + state.texture(a0 + 2, :) .* fr;
    amp = sqrt(10^(cl.powerDb/10) / sum(abs(tx).^2));        % per-sample power = powerDb
    SPEC = SPEC + fft(amp * z .* tex, nfft, 2) .* TX;
end

% --- targets ----------------------------------------------------------------------------
lossLong = sum(abs(tx).^2) * P.fft.nPulses;                     % coherent gain, upper bound
for k = 1:numel(S.targets)
    tg  = S.targets(k);
    G   = rsp_antenna_pattern(g.azDeg - tg.azDeg, P.antenna);   % one-way power
    a   = sqrt(10^(tg.powerDb/10)) * G;                        % two-way amplitude
    if max(a)^2 * lossLong < 0.01 * N0
        continue                                                % far below the noise
    end
    Rm  = tg.rangeM - tg.velocityMps * g.timeS;                 % range at each pulse
    tau = 2 * (Rm - P.radar.rangeOffsetM) / c * 1e6;            % delay [us]
    if any(tau * P.fs + txLen > nfft - 1) || any(Rm < 0)
        continue                                                % outside the record
    end
    ph  = exp(-1j*2*pi*mod(P.radar.fcMHz * tau, 1));            % carrier phase (Doppler)
    SPEC = SPEC + (a .* ph) .* exp(-1j*2*pi*tau * f) .* TX;
end

% --- receiver ------------------------------------------------------------------------------
rx    = ifft(SPEC, [], 2);
noise = complex(Z(state.cells + (1:R), :), Z(nZ/2 + state.cells + (1:R), :)).';
video = rx(:, 1:R) + sqrt(N0/2) * noise;
if S.blankTx
    video(:, 1:min(txLen, R)) = 0;
end
video = single(video);

if nargout >= 3
    truth = rsp_truth(P, S);
end
end
