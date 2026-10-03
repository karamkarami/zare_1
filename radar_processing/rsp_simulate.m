function [video, truth] = rsp_simulate(P, S)
%RSP_SIMULATE Synthetic baseband video for testing the chain.
%
%   [video, truth] = rsp_simulate(P, S)
%
%   P : chain parameters (P.fs and P.pulse define the transmitted signal)
%   S : scenario
%       .nPulses     number of pulses                       (350)
%       .nRange      range samples per pulse                (5469)
%       .noisePower  receiver noise power per sample        (1)
%       .targets     struct array: rangeCell (1-based), fdNorm (Doppler / PRF,
%                    -0.5..0.5), snrDb (per sample, before compression)
%       .clutter     struct: cnrDb (per sample), cells [first last],
%                    spreadNorm (Doppler std / PRF), [] = no clutter
%       .blankTx     true = zero the samples received during transmission
%       .seed        random seed ([] = do not reset)
%
%   video : nPulses x nRange single complex
%   truth : echo of S plus the transmitted signal and the blanked length

M = S.nPulses;
R = S.nRange;
if isfield(S, 'seed') && ~isempty(S.seed)
    rng(S.seed);
end

% --- transmitted signal: every pulse at its own delay ------------------------
nP = numel(P.pulse);
d  = zeros(1, nP);
s  = cell(1, nP);
for k = 1:nP
    s{k} = rsp_waveform(P.pulse(k), P.fs);
    d(k) = round(P.pulse(k).delayUs * P.fs);
end
txLen = max(d + cellfun(@numel, s));
tx    = zeros(txLen, 1);
for k = 1:nP
    tx(d(k) + (1:numel(s{k}))) = tx(d(k) + (1:numel(s{k}))) + s{k};
end
txEnergy = sum(abs(tx).^2);

% --- reflectivity (pulses x range cells) -------------------------------------
refl = zeros(M, R);
m    = (0:M-1)';
if isfield(S, 'clutter') && ~isempty(S.clutter)
    c   = S.clutter;
    rc  = max(1, c.cells(1)):min(R, c.cells(end));
    f   = [0:floor(M/2), -ceil(M/2)+1:-1]' / M;              % normalised Doppler
    psd = exp(-f.^2 / (2*c.spreadNorm^2));
    w   = (randn(M, numel(rc)) + 1j*randn(M, numel(rc))) / sqrt(2);
    cl  = ifft(bsxfun(@times, fft(w), sqrt(psd / mean(psd))));
    % per-sample CNR after the convolution with the transmitted signal
    refl(:, rc) = cl * sqrt(10^(c.cnrDb/10) * S.noisePower / txEnergy);
end
for k = 1:numel(S.targets)
    tg = S.targets(k);
    a  = sqrt(10^(tg.snrDb/10) * S.noisePower);
    refl(:, tg.rangeCell) = refl(:, tg.rangeCell) + a * exp(1j*2*pi*tg.fdNorm*m);
end

% --- echo = reflectivity (*) transmitted signal, along fast time ------------
nfft  = 2^nextpow2(R + txLen - 1);
rx    = ifft(fft(refl, nfft, 2) .* fft(tx, nfft).', [], 2);
video = rx(:, 1:R);
video = video + sqrt(S.noisePower/2) * (randn(M, R) + 1j*randn(M, R));
if S.blankTx
    video(:, 1:min(txLen, R)) = 0;
end
video = single(video);

truth          = S;
truth.tx       = tx;
truth.txLength = txLen;
end
