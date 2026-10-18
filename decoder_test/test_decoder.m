%% Simple decoder test for one PRI (short pulse + long pulse, two targets)
%
% Timing of one PRI:
%   | short Tx | gap1 | short listen | gap2 | long Tx | gap3 | long listen |
%
% Each pulse is an LFM (matched filter with a Taylor window) or a phase code
% (your mismatch filter, or a least-squares one, or the matched filter).
% Target 1 is at short range (seen in the short listening window) and
% target 2 is at long range (seen in the long listening window).

clear; clc; close all;

%% ---------------- Parameters ----------------
fs = 2e6;              % sampling frequency [Hz]
c  = 3e8;              % speed of light [m/s]
fc = 3e9;              % carrier frequency, only used for Doppler [Hz]

% timing [s]
gap1        = 2e-6;    % short Tx -> short listen
listenShort = 100e-6;  % short listening window
gap2        = 5e-6;    % short listen -> long Tx
gap3        = 2e-6;    % long Tx -> long listen
listenLong  = 800e-6;  % long listening window

% pulses: P(1) = short pulse, P(2) = long pulse
%  type 'code' : code    = phase code coefficients, one per chip (+1/-1 or exp(1j*phi))
%                chip    = chip width [s] (fs*chip samples per chip)
%                decoder = 'mmf' (mismatch filter) or 'matched'
%                mmf     = your mismatch filter taps ([] = least-squares design)
%                mmfRate = 'chip' (one tap per chip) or 'sample' (one tap per sample)
%                mmfLen  = mismatch length / code length, only for the LS design
%                          (3 -> 3*13 = 39 chips = 78 samples at 2 samples per chip)
%  type 'lfm'  : width [s], bw [Hz], window 'taylor' or 'none', nbar, sll [dB]
P(1).type    = 'code';
P(1).code    = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];   % Barker 13
P(1).chip    = 1e-6;
P(1).decoder = 'mmf';
P(1).mmf     = [];
P(1).mmfRate = 'chip';
P(1).mmfLen  = 3;

P(2).type    = 'lfm';
P(2).width   = 100e-6;
P(2).bw      = 1e6;
P(2).window  = 'taylor';
P(2).nbar    = 4;
P(2).sll     = -35;

% targets: range [m], radial velocity [m/s], SNR per sample before the decoder [dB]
tgt.range = [8e3   60e3];
tgt.vel   = [50   -100];
tgt.snr   = [5     -5];

thrDb = 13;            % detection threshold above the noise floor [dB]

%% ---------------- Pulses and decoders ----------------
names = {'short', 'long'};
for p = 1:2
    switch P(p).type
        case 'code'
            code = P(p).code(:);
            spc  = round(P(p).chip * fs);                 % samples per chip
            s    = kron(code, ones(spc, 1));              % transmitted pulse
            if strcmp(P(p).decoder, 'matched')
                h = conj(flipud(s));
            else
                hc = P(p).mmf(:);
                if isempty(hc)                            % least-squares mismatch filter
                    N = numel(code);
                    M = P(p).mmfLen * N;
                    A = zeros(N + M - 1, M);              % convolution matrix of the code
                    for k = 1:M
                        A(k:k+N-1, k) = code;
                    end
                    d = zeros(N + M - 1, 1);
                    d(ceil((N + M - 1) / 2)) = 1;         % wanted output: one peak in the middle
                    hc = A \ d;
                    P(p).mmf = hc;
                end
                if strcmp(P(p).mmfRate, 'chip')
                    h = kron(hc, ones(spc, 1));           % chip-rate taps held over the chip
                else
                    h = hc;
                end
            end
        case 'lfm'
            n = round(P(p).width * fs);
            t = ((0:n-1)' - (n-1)/2) / fs;
            s = exp(1j*pi*(P(p).bw/P(p).width) * t.^2);
            w = ones(n, 1);
            if strcmp(P(p).window, 'taylor')              % Taylor window (same as taylorwin)
                nb  = P(p).nbar;
                A   = acosh(10^(-P(p).sll/20)) / pi;
                sp2 = nb^2 / (A^2 + (nb - 0.5)^2);
                xi  = ((0:n-1)' - (n-1)/2) / n;
                for m = 1:nb-1
                    num = prod(1 - (m^2/sp2) ./ (A^2 + ((1:nb-1) - 0.5).^2));
                    den = prod(1 - m^2 ./ [1:m-1, m+1:nb-1].^2);
                    w   = w + (-1)^(m+1) * num / den * cos(2*pi*m*xi);
                end
            end
            h = conj(flipud(s)) .* w;
    end
    h = h / norm(h);                                      % unit noise gain: output in dB = SNR

    r = conv(s, h);                                       % noise-free compressed pulse
    [pk, k] = max(abs(r));
    lo = k; while lo > 1         && abs(r(lo-1)) < abs(r(lo)), lo = lo - 1; end
    hi = k; while hi < numel(r)  && abs(r(hi+1)) < abs(r(hi)), hi = hi + 1; end
    side = abs(r([1:lo-1, hi+1:end]));

    pulse(p).s    = s;
    pulse(p).h    = h;
    pulse(p).r    = r;
    pulse(p).lag  = k - 1;                                % peak lag for an echo starting at 0
    pulse(p).gain = 20*log10(pk);                         % SNR gain of the decoder [dB]
    pulse(p).loss = 10*log10(numel(s)) - pulse(p).gain;   % loss against a matched filter [dB]
    pulse(p).psl  = 20*log10(max([side; eps]) / pk);      % peak sidelobe level [dB]
    fprintf('%-5s pulse: %s, %d samples, decoder %d taps, gain %.1f dB, loss %.2f dB, PSL %.1f dB\n', ...
        names{p}, P(p).type, numel(s), numel(h), pulse(p).gain, pulse(p).loss, pulse(p).psl);
end

%% ---------------- Timeline of one PRI ----------------
iTx(1) = 0;                                               % start samples (0 = PRI start)
iRx(1) = iTx(1) + numel(pulse(1).s) + round(gap1 * fs);
nRx(1) = round(listenShort * fs);
iTx(2) = iRx(1) + nRx(1) + round(gap2 * fs);
iRx(2) = iTx(2) + numel(pulse(2).s) + round(gap3 * fs);
nRx(2) = round(listenLong * fs);
nPRI   = iRx(2) + nRx(2);

tx = zeros(nPRI, 1);
tx(iTx(1) + (1:numel(pulse(1).s))) = pulse(1).s;
tx(iTx(2) + (1:numel(pulse(2).s))) = pulse(2).s;

%% ---------------- Received signal ----------------
t  = (0:nPRI-1)' / fs;
rx = (randn(nPRI, 1) + 1j*randn(nPRI, 1)) / sqrt(2);     % noise, power 1
for k = 1:numel(tgt.range)
    d    = round(2 * tgt.range(k) / c * fs);             % delay [samples]
    fd   = 2 * tgt.vel(k) * fc / c;                      % Doppler [Hz]
    a    = 10^(tgt.snr(k)/20) * exp(1j*2*pi*rand);
    echo = [zeros(d, 1); tx(1:nPRI-d)];
    rx   = rx + a * echo .* exp(1j*2*pi*fd*t);
end

%% ---------------- Decoder in each listening window ----------------
figure('Name', 'Decoder test');
for p = 1:2
    seg = rx(iRx(p) + (1:nRx(p)));                        % samples of this listening window
    y   = conv(seg, pulse(p).h);
    y   = y(pulse(p).lag + (1:nRx(p)));                   % sample i = echo starting at sample i
    rng = c/2 * ((iRx(p) - iTx(p)) + (0:nRx(p)-1)') / fs; % range of each output sample [m]

    ydb     = 20*log10(abs(y) + eps);
    noiseDb = 10*log10(median(abs(y).^2) / log(2));       % noise floor (robust to targets)
    thr     = noiseDb + thrDb;
    isPeak  = [false; ydb(2:end-1) > ydb(1:end-2) & ydb(2:end-1) >= ydb(3:end); false];
    det     = find(isPeak & ydb > thr);

    fprintf('\n%s window: %.2f - %.2f km, noise floor %.1f dB\n', names{p}, rng(1)/1e3, rng(end)/1e3, noiseDb);
    for k = 1:numel(tgt.range)
        if tgt.range(k) >= rng(1) && tgt.range(k) <= rng(end)
            fprintf('   target %d: range %.3f km, expected SNR %.1f dB\n', ...
                k, tgt.range(k)/1e3, tgt.snr(k) + pulse(p).gain);
        end
    end
    for k = det'
        fprintf('   detection: range %.3f km, SNR %.1f dB\n', rng(k)/1e3, ydb(k) - noiseDb);
    end
    if isempty(det)
        fprintf('   no detection\n');
    end

    % noise-free compressed pulse
    subplot(2, 2, p);
    rdb  = 20*log10(abs(pulse(p).r) / max(abs(pulse(p).r)) + eps);
    lags = ((0:numel(rdb)-1)' - pulse(p).lag) / fs * 1e6;
    plot(lags, rdb); grid on; ylim([-80 5]);
    xlabel('lag [\mus]'); ylabel('dB');
    title({sprintf('%s pulse (%s)', names{p}, P(p).type), ...
           sprintf('PSL %.1f dB, loss %.2f dB', pulse(p).psl, pulse(p).loss)});

    % decoder output of the listening window
    subplot(2, 2, p + 2);
    plot(rng/1e3, ydb - noiseDb); hold on; grid on;
    plot(rng([1 end])/1e3, [thrDb thrDb], 'r--');
    plot(rng(det)/1e3, ydb(det) - noiseDb, 'ko');
    for k = 1:numel(tgt.range)
        if tgt.range(k) >= rng(1) && tgt.range(k) <= rng(end)
            plot(tgt.range(k)/1e3 * [1 1], [-10 60], 'g:');
        end
    end
    xlabel('range [km]'); ylabel('SNR [dB]'); ylim([-10 max(40, max(ydb - noiseDb) + 5)]);
    title(sprintf('%s listening window', names{p}));
end

% one PRI: |rx| with the transmit and listening windows
figure('Name', 'One PRI');
plot(t*1e6, 20*log10(abs(rx) + eps)); hold on; grid on;
yl = [-30 30];
edges = [iTx(1), iTx(1)+numel(pulse(1).s); iRx(1), iRx(1)+nRx(1); ...
         iTx(2), iTx(2)+numel(pulse(2).s); iRx(2), iRx(2)+nRx(2)] / fs * 1e6;
cols  = {'r', 'g', 'r', 'g'};
for k = 1:4
    plot(edges(k, [1 1 2 2]), yl([1 2 2 1]), cols{k});
end
ylim(yl); xlabel('time [\mus]'); ylabel('|rx| [dB]');
title('One PRI: red = Tx, green = listening windows');
