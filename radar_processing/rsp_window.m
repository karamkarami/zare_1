function w = rsp_window(spec, N, param)
%RSP_WINDOW Amplitude window of length N (column vector, peak = 1).
%
%   w = rsp_window(spec, N, param)
%
%   spec  : 'none' | 'rect' | 'hamming' | 'hann' | 'blackman' |
%           'blackmanharris' | 'kaiser' | 'taylor'
%           or a numeric vector of length N, or a function handle @(N)
%   param : kaiser -> beta (default 6)
%           taylor -> [nbar sllDb] (default [4 -35])
%
%   No toolbox is required.

if nargin < 3
    param = [];
end
N = double(N);
if N <= 0
    w = zeros(0, 1);
    return
end

if isnumeric(spec)
    w = double(spec(:));
    if numel(w) ~= N
        error('rsp_window:size', 'Numeric window has %d taps, expected %d.', numel(w), N);
    end
    return
end
if isa(spec, 'function_handle')
    w = double(spec(N));
    w = w(:);
    return
end
if N == 1
    w = 1;
    return
end

n = (0:N-1)';
x = 2*pi*n/(N-1);
switch lower(spec)
    case {'none', 'rect', 'rectangular', 'boxcar'}
        w = ones(N, 1);
    case 'hamming'
        w = 0.54 - 0.46*cos(x);
    case {'hann', 'hanning'}
        w = 0.5 - 0.5*cos(x);
    case 'blackman'
        w = 0.42 - 0.5*cos(x) + 0.08*cos(2*x);
    case 'blackmanharris'
        w = 0.35875 - 0.48829*cos(x) + 0.14128*cos(2*x) - 0.01168*cos(3*x);
    case 'kaiser'
        if isempty(param), param = 6; end
        r = 2*n/(N-1) - 1;
        w = besseli(0, param*sqrt(max(0, 1 - r.^2))) / besseli(0, param);
    case 'taylor'
        if isempty(param), param = [4 -35]; end
        w = taylorWindow(N, param(1), param(2));
    otherwise
        error('rsp_window:type', 'Unknown window "%s".', spec);
end
w = w / max(w);
end

function w = taylorWindow(N, nbar, sll)
% Taylor window, same definition as taylorwin (Signal Processing Toolbox).
A   = acosh(10^(-sll/20)) / pi;
sp2 = nbar^2 / (A^2 + (nbar - 0.5)^2);
m   = 1:nbar-1;
Fm  = zeros(size(m));
for k = m
    num    = prod(1 - (k^2/sp2) ./ (A^2 + (m - 0.5).^2));
    others = m(m ~= k);
    den    = prod(1 - k^2 ./ others.^2);
    Fm(k)  = ((-1)^(k+1) / 2) * num / den;
end
x = ((0:N-1)' - (N-1)/2) / N;
w = 1 + 2*cos(2*pi*x*m) * Fm';
end
