function [video, info, P, truth] = mc_synthetic(P, S)
%MC_SYNTHETIC Test video of the two-pulse radar: noise, clutter and moving targets.
%
%   [video, info, P, truth] = mc_synthetic(P)      default scene
%   [video, info, P, truth] = mc_synthetic(P, S)   your scene (fields below)
%
%   Uses the pulses of P (lengths, delays, frequencies). A pulse without a code
%   gets a random binary code; a pulse without decoder coefficients gets a
%   least-squares mismatched filter of coefLength taps (unit noise gain, so both
%   pulses have the same noise floor). Both are written back into P, so the
%   chain decodes the test signal with them.
%
%   S.nPulses, S.nCells   size of the video                      (512 x 2000)
%   S.noiseDb             noise power per sample                 (0 dB)
%   S.clutterDb           zero-Doppler clutter power per sample  (20 dB)
%   S.clutterCells        clutter up to this cell                (900)
%   S.clutterSpread       pulse-to-pulse fluctuation of the clutter (1e-3)
%   S.targets             struct array: cell, bin (Doppler bin 1..nfft), powerDb
%   S.blank               receiver off while transmitting        (true)
%   S.seed                random seed                            (1)
%
%   Echo model: range cell r reflects every pulse k; its echo starts at sample
%   r + txDelay_k, as mc_filter assumes.

if nargin < 2
    S = struct();
end
S = setDefault(S, 'nPulses', 512);
S = setDefault(S, 'nCells', 2000);
S = setDefault(S, 'noiseDb', 0);
S = setDefault(S, 'clutterDb', 20);
S = setDefault(S, 'clutterCells', 900);
S = setDefault(S, 'clutterSpread', 1e-3);
S = setDefault(S, 'targets', struct('cell', {300, 700, 1200}, 'bin', {5, 9, 12}, ...
                                    'powerDb', {0, -5, -10}));
S = setDefault(S, 'blank', true);
S = setDefault(S, 'seed', 1);
rng(S.seed);

% ---- transmitted frame: every pulse at its delay and frequency --------------------------------
fs = P.fsMHz;
tx = [];
for k = 1:numel(P.pulse)
    pk  = P.pulse(k);
    nCh = round(pk.widthUs * pk.bwMHz);
    if isempty(pk.code)
        pk.code = sign(randn(nCh, 1));
        pk.code(pk.code == 0) = 1;
        P.pulse(k).code = pk.code;
    end
    spc = fs / pk.bwMHz;
    s   = pk.code(min(floor((0:round(numel(pk.code) * spc) - 1)' / spc) + 1, numel(pk.code)));
    if isempty(pk.coef) && isempty(pk.coefFile)
        h = mc_mismatch(pk.code, pk.coefLength, spc);
        h = h / norm(h);
        if strcmpi(P.decoder.coefForm, 'correlate')
            h = conj(flipud(h));
        end
        P.pulse(k).coef = h;
    end
    d = round(pk.txDelayUs * fs);
    n = (0:numel(s) - 1)';
    if numel(tx) < d + numel(s)
        tx(d + numel(s), 1) = 0;
    end
    tx(d + (1:numel(s))) = tx(d + (1:numel(s))) + s .* exp(1j*2*pi*pk.freqMHz/fs*(d + n));
end

% ---- reflectivity: pulses x range cells -------------------------------------------------------
nP = S.nPulses;
nC = S.nCells;
p  = (0:nP - 1)';
refl = zeros(nP, nC);
c  = 1:min(S.clutterCells, nC);
g  = sqrt(10^(S.clutterDb / 10) / 2) * (randn(1, numel(c)) + 1j*randn(1, numel(c)));
refl(:, c) = repmat(g, nP, 1) .* (1 + S.clutterSpread * (randn(nP, numel(c)) + 1j*randn(nP, numel(c))));
truth = S.targets;
for i = 1:numel(truth)
    fd = (truth(i).bin - 1) / P.doppler.nfft;                  % cycles per pulse
    refl(:, truth(i).cell) = refl(:, truth(i).cell) + ...
        sqrt(10^(truth(i).powerDb / 10)) * exp(1j*2*pi*fd*p);
end

% ---- video = reflectivity (*) transmitted frame, + noise ----------------------------------------
video = filter(tx, 1, refl, [], 2);                            % echo of cell r starts at sample r
video = video + sqrt(10^(S.noiseDb / 10) / 2) * (randn(nP, nC) + 1j*randn(nP, nC));
if S.blank
    video(:, 1:min(numel(tx), nC)) = 0;
end
video = single(video);

info.pulseSeq   = 1000 + (1:nP)';
info.azimuthDeg = (0:nP - 1)' * 360 / 3560;
end

function S = setDefault(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end
