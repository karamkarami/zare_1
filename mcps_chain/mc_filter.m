function f = mc_filter(pulse, P)
%MC_FILTER Decoder filter of one pulse, its range alignment and quality figures.
%
%   f = mc_filter(P.pulse(k), P)
%
%   f.h            filter impulse response at fs (column): z = conv(x, h)
%   f.nSamples     pulse length [samples]
%   f.nCoef        number of decoder coefficients
%   f.lag          lag of the main peak for an echo that starts at sample 0
%   f.shift        lag + transmit delay [samples]: output cell r = z(r + shift)
%   f.lagFrom      'code' (peak measured with pulse.code), 'centred' (assumes a
%                  mismatched filter centred on the code) or 'user' (pulse.alignLag)
%   f.placeholder  true when no coefficients were given
%   f.pslDb, f.lossDb   peak sidelobe, and gain against a matched filter (<= 0 dB,
%                  -3 = 3 dB mismatch loss); both need pulse.code

fs = P.fsMHz;
N  = round(pulse.widthUs * fs);                       % pulse length in samples
s  = codeSamples(pulse, fs, N);                       % [] when no code is given

% ---- coefficients --------------------------------------------------------------------
coef = pulse.coef(:);
if isempty(coef) && ~isempty(pulse.coefFile)
    coef = loadCoef(pulse.coefFile);
end
placeholder = isempty(coef);
form = lower(P.decoder.coefForm);
if placeholder
    % Placeholder: matched filter of the code (or a plain pulse of N samples),
    % centred in coefLength taps like a mismatched filter. NOT your decoder.
    if isempty(s)
        g = ones(N, 1);
    else
        g = conj(flipud(s));
    end
    L    = max(pulse.coefLength, N);
    a    = floor((L - N) / 2);
    coef = [zeros(a, 1); g; zeros(L - N - a, 1)];
    form = 'fir';
elseif ~isempty(pulse.coefLength) && numel(coef) ~= pulse.coefLength
    warning('mc_filter:length', 'Pulse "%s": %d coefficients, %d expected.', ...
            pulse.name, numel(coef), pulse.coefLength);
end

switch form
    case 'fir',       h = coef;
    case 'correlate', h = conj(flipud(coef));
    otherwise, error('mc_filter:form', 'Unknown P.decoder.coefForm "%s".', form);
end
h = double(h);
if strcmpi(P.decoder.normalize, 'noise')
    h = h / norm(h);
end
h = h * pulse.gain;

% ---- optional low-pass (frequency channel of this pulse) -------------------------------
lpDelay = 0;
if P.decoder.bandFilter
    lp = mc_lowpass(P.decoder.bandFactor * pulse.bwMHz / 2, fs, P.decoder.bandTaps);
    h  = conv(lp, h);
    lpDelay = (numel(lp) - 1) / 2;
end

% ---- range alignment ---------------------------------------------------------------------
f.pslDb  = NaN;
f.lossDb = NaN;
if ~isempty(pulse.alignLag)
    lag  = pulse.alignLag;
    from = 'user';
elseif ~isempty(s)
    r = conv(s, h);                                   % compressed pulse, echo at sample 0
    [~, i] = max(abs(r));
    lag  = i - 1;
    from = 'code';
else
    lag  = round((N + numel(coef)) / 2 - 1 + lpDelay);  % centred mismatched filter
    from = 'centred';
end
if ~isempty(s)
    [f.pslDb, f.lossDb] = quality(conv(s, h), s, h);
end

f.h           = h;
f.nSamples    = N;
f.nCoef       = numel(coef);
f.lag         = lag;
f.shift       = lag + round(pulse.txDelayUs * fs);
f.lagFrom     = from;
f.placeholder = placeholder;
if placeholder
    warning('mc_filter:placeholder', ['Pulse "%s" has no decoder coefficients: a ' ...
            'placeholder filter is used. Put yours in mc_params.m.'], pulse.name);
end
end

function s = codeSamples(pulse, fs, N)
% chips -> samples at fs (each chip held for fs/bw samples)
s = [];
if isempty(pulse.code)
    return
end
c   = double(pulse.code(:));
spc = fs / pulse.bwMHz;                               % samples per chip
k   = floor((0:round(numel(c) * spc) - 1)' / spc) + 1;
s   = c(min(k, numel(c)));
if numel(s) ~= N
    warning('mc_filter:code', 'Pulse "%s": code gives %d samples, widthUs gives %d.', ...
            pulse.name, numel(s), N);
end
end

function c = loadCoef(file)
% .mat: first numeric variable; text: one value per line or two columns [real imag]
[~, ~, ext] = fileparts(file);
if strcmpi(ext, '.mat')
    S = load(file);
    names = fieldnames(S);
    c = [];
    for i = 1:numel(names)
        if isnumeric(S.(names{i}))
            c = S.(names{i});
            break
        end
    end
else
    c = load(file);
end
if isempty(c)
    error('mc_filter:file', 'No coefficients found in "%s".', file);
end
if size(c, 2) == 2 && isreal(c)
    c = complex(c(:, 1), c(:, 2));
end
c = c(:);
end

function [pslDb, lossDb] = quality(r, s, h)
% PSL: highest sample outside the main lobe (down to its first nulls)
a = abs(r);
[pk, k] = max(a);
lo = k;
while lo > 1 && a(lo - 1) < a(lo)
    lo = lo - 1;
end
hi = k;
while hi < numel(a) && a(hi + 1) < a(hi)
    hi = hi + 1;
end
side = a([1:lo-1, hi+1:end]);
pslDb  = 20*log10(max([side; eps]) / pk);
lossDb = 20*log10(pk / (norm(s) * norm(h)));          % 0 dB = matched filter
end
