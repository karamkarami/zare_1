function [Y, info] = rsp_matched_filter(X, P)
%RSP_MATCHED_FILTER Pulse compression of a multi-pulse (e.g. short + long) waveform.
%
%   [Y, info] = rsp_matched_filter(X, P)
%
%   X    : video, pulses x range samples (complex baseband)
%   P    : parameter struct (uses P.fs, P.pulse, P.mf)
%   Y    : decoder output, pulses x range cells (same size as X)
%   info : .replica       cell, matched-filter taps of each pulse
%          .delay         transmit delay of each pulse [samples]
%          .length        length of each pulse [samples]
%          .switchCell    number of near cells taken from the short pulse
%          .shortPulse    index of the short pulse
%          .longPulse     index of the long pulse
%          .nfft          FFT length used
%          .each          cell, MF output of each pulse (if P.mf.keepEach)
%
%   Range alignment: output cell r (1-based) holds the response of a
%   scatterer whose echo starts at sample r + delay of that pulse, i.e. all
%   pulses are aligned to the same range grid before they are combined.
%
%   The correlation is done in the frequency domain with one forward FFT of
%   the data and one inverse FFT per pulse (no circular wrap-around).

[M, R] = size(X);
cls    = class(X);
mf     = P.mf;
nP     = numel(P.pulse);

% --- replicas -----------------------------------------------------------
h = cell(1, nP);
L = zeros(1, nP);
d = zeros(1, nP);
for k = 1:nP
    pk = P.pulse(k);
    s  = rsp_waveform(pk, P.fs);
    w  = rsp_window(pk.window, numel(s), pk.windowParam);
    hk = s .* w;
    switch lower(mf.norm)
        case 'noise', hk = hk / norm(hk);              % noise power gain = 1
        case 'peak',  hk = hk / abs(sum(s .* conj(hk)));  % matched echo peak = |amplitude|
        case 'none'
        otherwise, error('rsp_matched_filter:norm', 'Unknown P.mf.norm "%s".', mf.norm);
    end
    if isfield(pk, 'gainDb') && ~isempty(pk.gainDb)
        hk = hk * 10^(pk.gainDb/20);
    end
    h{k} = hk;
    L(k) = numel(hk);
    d(k) = round(pk.delayUs * P.fs);
end

[~, iLong]  = max(L);
[~, iShort] = min(L);
if isempty(mf.switchCell)
    % Echoes of cells closer than this overlap the transmission of the
    % long pulse (eclipsing): use the short pulse there.
    switchCell = max(d + L) - d(iLong);
else
    switchCell = mf.switchCell;
end
switchCell = min(max(round(switchCell), 0), R);

info = struct('replica', {h}, 'delay', d, 'length', L, 'switchCell', switchCell, ...
              'shortPulse', iShort, 'longPulse', iLong, 'nfft', 0, 'each', {{}});

if ~mf.enable
    Y = X;
    return
end

% --- which pulses are needed ----------------------------------------------
combine = lower(mf.combine);
if nP == 1
    sel  = 1;
    need = 1;
elseif strcmp(combine, 'stitch')
    sel  = [];
    need = unique([iShort iLong]);
elseif strncmp(combine, 'pulse', 5)
    sel = str2double(combine(6:end));
    if isnan(sel) || sel < 1 || sel > nP
        error('rsp_matched_filter:combine', 'P.mf.combine = "%s" is not a valid pulse.', mf.combine);
    end
    need = sel;
else
    error('rsp_matched_filter:combine', 'Unknown P.mf.combine "%s".', mf.combine);
end
if mf.keepEach
    need = 1:nP;
end

% --- frequency-domain correlation ------------------------------------------
nfft      = fastLength(R + max(L) - 1);
info.nfft = nfft;
FX        = fft(X, nfft, 2);
each      = cell(1, nP);
cells     = 0:R-1;
for k = need
    Hk  = conj(fft(cast(h{k}, cls), nfft)).';           % 1 x nfft
    yk  = ifft(FX .* Hk, [], 2);                         % y[n] = sum x[n+m] h*[m]
    n   = cells + d(k) + mf.rangeOffset;                 % lag of each output cell
    ok  = n >= 0 & n <= R-1;
    out = complex(zeros(M, R, cls));
    out(:, ok) = yk(:, n(ok) + 1);
    each{k} = out;
end
clear FX yk

% --- combine ---------------------------------------------------------------
if isempty(sel)
    Y = each{iLong};
    Y(:, 1:switchCell) = each{iShort}(:, 1:switchCell);
else
    Y = each{sel};
end
if mf.keepEach
    info.each = each;
end
end

function n = fastLength(n)
% Smallest integer >= n whose prime factors are 2, 3 and 5 (fast FFT size).
while true
    m = n;
    for p = [2 3 5]
        while mod(m, p) == 0
            m = m / p;
        end
    end
    if m == 1
        return
    end
    n = n + 1;
end
end
