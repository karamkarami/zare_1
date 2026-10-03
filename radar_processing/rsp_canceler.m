function [Y, idx, w] = rsp_canceler(X, C, idx)
%RSP_CANCELER N-pulse MTI canceler along slow time (default: 3-pulse).
%
%   [Y, idx, w] = rsp_canceler(X, C, idx)
%
%   X   : pulses x range (complex)
%   C   : P.canceler (see rsp_default_params)
%   idx : pulse index of each row of X (default 1:size(X,1))
%   Y   : canceler output, y(m) = sum_k w(k) x(m-k+1)
%         3-pulse: y(m) = x(m) - 2 x(m-1) + x(m-2)
%   idx : pulse index of each row of Y
%   w   : slow-time taps used

if nargin < 3 || isempty(idx)
    idx = (1:size(X, 1))';
end
idx = idx(:);

w = rsp_canceler_taps(C);
if ~C.enable
    Y = X;
    return
end

nT = numel(w) - 1;
if size(X, 1) <= nT
    error('rsp_canceler:pulses', 'Need more than %d pulses for a %d-tap canceler.', nT, nT + 1);
end
Y = filter(cast(w, class(X)), 1, X, [], 1);

switch lower(C.output)
    case 'same'
        Y(1:nT, :) = 0;                 % filter transient
    case 'valid'
        Y   = Y(nT+1:end, :);
        idx = idx(nT+1:end);
    otherwise
        error('rsp_canceler:output', 'Unknown P.canceler.output "%s".', C.output);
end
end
