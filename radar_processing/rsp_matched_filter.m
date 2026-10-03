function [Y, info] = rsp_matched_filter(X, P)
%RSP_MATCHED_FILTER Pulse compression of a two-pulse (short + long) waveform.
%
%   [Y, info] = rsp_matched_filter(X, P)
%
%   X    : video, pulses x range samples (complex baseband)
%   P    : parameters (uses P.fs, P.pulse, P.mf)
%   Y    : decoder output, pulses x range cells (same size as X)
%   info : .dec          rsp_decoder output of each pulse (h, lag, PSL, loss ...)
%          .replica      cell, correlation reference of each pulse
%          .delay        transmit delay of each pulse [samples]
%          .length       length of each pulse [samples]
%          .switchCell   cells 1..switchCell come from the short pulse
%          .validCells   [first last] cells whose echo is received in full
%                        (transmit eclipsing at the start, record end at the end)
%          .shortPulse, .longPulse   pulse indices
%          .nfft         FFT length
%          .each         cell, output of each pulse (if P.mf.keepEach)
%
%   Each pulse is either an LFM or a phase code, decoded with the matched
%   filter or with your own decoder taps (see rsp_decoder). The correlation
%   runs in the frequency domain: one forward FFT of the data and one
%   inverse FFT per pulse, without circular wrap-around.
%
%   Range alignment: output cell r (1-based) holds the main peak of a
%   scatterer whose echo of pulse k starts at sample r-1 + delay_k, for every
%   pulse and every decoder (the decoder peak lag is removed).

[M, R] = size(X);
cls    = class(X);
mf     = P.mf;
nP     = numel(P.pulse);

dec = cell(1, nP);
h   = cell(1, nP);
L   = zeros(1, nP);
d   = zeros(1, nP);
for k = 1:nP
    dec{k} = rsp_decoder(P.pulse(k), P.fs, mf.norm, mf);
    h{k}   = dec{k}.h;
    L(k)   = numel(dec{k}.tx);
    d(k)   = round(P.pulse(k).delayUs * P.fs);
end

[~, iLong]  = max(L);
[~, iShort] = min(L);
if isempty(mf.switchCell)
    % Echoes from cells closer than this overlap the transmission of the
    % long pulse (eclipsing): use the short pulse there.
    switchCell = max(d + L) - d(iLong);
else
    switchCell = mf.switchCell;
end
switchCell = min(max(round(switchCell), 0), R);

% cells received in full: the short pulse echo clears the transmission
% at the start, the long pulse echo still fits in the record at the end
iFirst = iShort;
iLast  = iLong;
if nP > 1 && strncmpi(mf.combine, 'pulse', 5)
    iFirst = str2double(mf.combine(6:end));
    iLast  = iFirst;
end
first = max(d + L) - d(iFirst) + 1;
last  = R - L(iLast) - d(iLast) + 1;
validCells = [min(max(first, 1), R) min(max(last, 1), R)];

info = struct('dec', {dec}, 'replica', {h}, 'delay', d, 'length', L, ...
              'switchCell', switchCell, 'validCells', validCells, ...
              'shortPulse', iShort, 'longPulse', iLong, 'nfft', 0, 'each', {{}});
if ~mf.enable
    Y = X;
    return
end

% --- pulses to compute -----------------------------------------------------------
combine = lower(mf.combine);
if nP == 1
    sel = 1;
elseif strcmp(combine, 'stitch')
    sel = [];
elseif strncmp(combine, 'pulse', 5)
    sel = str2double(combine(6:end));
    if isnan(sel) || sel < 1 || sel > nP
        error('rsp_matched_filter:combine', 'P.mf.combine = "%s" is not a valid pulse.', mf.combine);
    end
else
    error('rsp_matched_filter:combine', 'Unknown P.mf.combine "%s".', mf.combine);
end
if mf.keepEach
    need = 1:nP;
elseif isempty(sel)
    need = unique([iShort iLong]);
else
    need = sel;
end

% --- frequency-domain correlation ---------------------------------------------------
Lh        = max(cellfun(@numel, h));
nfft      = fastLength(R + Lh - 1);
info.nfft = nfft;
FX        = fft(X, nfft, 2);
each      = cell(1, nP);
cells     = 0:R-1;
for k = need
    Hk  = conj(fft(cast(h{k}, cls), nfft)).';           % 1 x nfft
    yk  = ifft(FX .* Hk, [], 2);                         % y(n) = sum x(n+m) h*(m)
    n   = cells + d(k) + dec{k}.lag + mf.rangeOffset;    % lag read by each cell
    ok  = n >= -(numel(h{k}) - 1) & n <= R - 1;
    out = complex(zeros(M, R, cls));
    out(:, ok) = yk(:, mod(n(ok), nfft) + 1);
    each{k} = out;
end
clear FX yk

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
