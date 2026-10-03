function dec = rsp_decoder(pulse, fs, norm)
%RSP_DECODER Receive filter (decoder) of one pulse and its quality figures.
%
%   dec = rsp_decoder(pulse, fs, norm)
%
%   pulse : one entry of P.pulse
%   fs    : sampling frequency [MHz]
%   norm  : 'noise' (unit noise gain) | 'peak' (unit gain for a matched echo) | 'none'
%
%   dec.tx      transmitted samples (column)
%   dec.h       correlation reference: y(n) = sum_m x(n+m) conj(h(m))
%   dec.lag     lag of the main peak for an echo starting at lag 0 [samples];
%               the matched filter shifts its output by this lag
%   dec.lossDb  mismatch loss against a matched filter (0 dB = matched)
%   dec.pslDb   peak sidelobe level of the compressed pulse
%   dec.islDb   integrated sidelobe level
%   dec.response compressed pulse (noise free), dec.lags its lags
%
%   Decoder
%     pulse.decoder = []          matched filter (replica, optionally windowed)
%     pulse.decoder = vector      your decoder coefficients:
%        pulse.decoderForm 'fir'       FIR taps applied by convolution
%                                      (matched = conj(fliplr(code)))
%                          'reference' correlation reference (matched = code)
%        pulse.decoderRate 'chip'      one tap per chip (held over the chip)
%                          'sample'    one tap per sample at fs
%     pulse.decoderLag = []       automatic alignment on the main peak, or a
%                                 lag in samples
%     pulse.window                amplitude weighting of h (LFM sidelobes)

if nargin < 3 || isempty(norm)
    norm = 'noise';
end
[tx, spc] = rsp_waveform(pulse, fs);

d = getField(pulse, 'decoder', []);
if isempty(d)
    h = tx;                                            % matched filter
else
    d = double(d(:));
    if strcmpi(getField(pulse, 'decoderRate', 'chip'), 'chip') && ...
            strcmpi(getField(pulse, 'type', 'lfm'), 'code') && ...
            ~strcmpi(getField(pulse, 'codeRate', 'chip'), 'sample')
        L    = round(numel(d) * spc);
        tap  = floor((0:L-1)' / spc) + 1;
        d    = d(min(tap, numel(d)));                   % chip-rate taps held over the chip
    end
    fc = getField(pulse, 'fcMHz', 0);
    if fc ~= 0                                           % decoder follows the pulse offset
        n = (0:numel(d)-1)';
        if strcmpi(getField(pulse, 'decoderForm', 'fir'), 'fir')
            n = n - (numel(d) - 1);                      % taps run backwards in time
        end
        d = d .* exp(1j*2*pi*fc*n/fs);
    end
    switch lower(getField(pulse, 'decoderForm', 'fir'))
        case 'fir',       h = conj(flipud(d));
        case 'reference', h = d;
        otherwise
            error('rsp_decoder:form', 'Unknown decoderForm "%s".', pulse.decoderForm);
    end
end
h = h .* rsp_window(getField(pulse, 'window', 'none'), numel(h), ...
                    getField(pulse, 'windowParam', []));

% noise-free compressed pulse: lags -(Lh-1) .. Ls-1
y    = conv(tx, conj(flipud(h)));
lags = (0:numel(y)-1)' - (numel(h) - 1);
[pk, k] = max(abs(y));
lag  = getField(pulse, 'decoderLag', []);
if isempty(lag)
    lag = lags(k);
end

lossDb = 20*log10(pk / (norm2(tx) * norm2(h)));       % before any scaling

switch lower(norm)
    case 'noise', g = 1 / norm2(h);
    case 'peak',  g = 1 / pk;
    case 'none',  g = 1;
    otherwise, error('rsp_decoder:norm', 'Unknown norm "%s".', norm);
end
g = g * 10^(getField(pulse, 'gainDb', 0)/20);
h = h * conj(g);                                       % y scales by g
y = y * g;

% quality: sidelobes = everything outside the first nulls around the peak
a  = abs(y);
lo = k;
while lo > 1 && a(lo - 1) < a(lo)
    lo = lo - 1;
end
hi = k;
while hi < numel(a) && a(hi + 1) < a(hi)
    hi = hi + 1;
end
main = false(size(a));
main(lo:hi) = true;
side  = a(~main);
pkAbs = abs(y(k));
dec.tx       = tx;
dec.h        = h;
dec.lag      = lag;
dec.spc      = spc;
dec.lossDb   = lossDb;
dec.pslDb    = 20*log10(max([side; eps]) / pkAbs);
dec.islDb    = 10*log10(sum(side.^2) / pkAbs^2 + eps);
dec.response = y;
dec.lags     = lags;
end

function n = norm2(x)
n = sqrt(sum(abs(x(:)).^2));
end

function v = getField(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
