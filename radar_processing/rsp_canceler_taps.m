function w = rsp_canceler_taps(C)
%RSP_CANCELER_TAPS Slow-time taps of the MTI canceler (row vector).
%
%   w = rsp_canceler_taps(C)       C = P.canceler
%
%   Binomial taps with alternating sign for C.order pulses
%   (2-pulse [1 -1], 3-pulse [1 -2 1], 4-pulse [1 -3 3 -1] ...), or
%   C.weights when given. With C.normalize the taps have unit norm.

if isfield(C, 'weights') && ~isempty(C.weights)
    w = double(C.weights(:)).';
else
    k = 0:C.order-1;
    w = (-1).^k .* arrayfun(@(i) nchoosek(C.order - 1, i), k);
end
if C.normalize
    w = w / norm(w);
end
end
