function w = mc_window(name, N)
%MC_WINDOW Window as a column vector (no toolbox needed).
%
%   w = mc_window(name, N)    name: 'rect' | 'hann' | 'hamming' | 'blackman',
%                             or a numeric vector of N weights

if isnumeric(name)
    w = double(name(:));
    if numel(w) ~= N
        error('mc_window:length', 'Window has %d values, %d expected.', numel(w), N);
    end
    return
end
if N == 1
    w = 1;
    return
end
x = 2 * pi * (0:N-1)' / (N - 1);
switch lower(name)
    case {'rect', 'none'}, w = ones(N, 1);
    case 'hann',           w = 0.5  - 0.5  * cos(x);
    case 'hamming',        w = 0.54 - 0.46 * cos(x);
    case 'blackman',       w = 0.42 - 0.5  * cos(x) + 0.08 * cos(2 * x);
    otherwise, error('mc_window:name', 'Unknown window "%s".', name);
end
end
