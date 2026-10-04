function [Y, rows] = mc_canceler(X, P, rows)
%MC_CANCELER MTI canceler along slow time (from pulse to pulse).
%
%   [Y, rows] = mc_canceler(X, P, rows)
%
%   X : pulses x range cells (decoder output)
%   Y : same size. With P.canceler.coef = [1 -2 1] (3-pulse canceler):
%           Y(n, r) = X(n, r) - 2 X(n-1, r) + X(n-2, r)
%       The first numel(coef)-1 rows are 0 (not enough pulses) and are
%       marked rows.valid = false.

if ~P.canceler.enable
    Y = X;
    return
end
b = P.canceler.coef(:);
K = numel(b);
Y = complex(zeros(size(X), class(X)));
for k = 1:K                                   % b(k) weights the pulse k-1 back
    Y(K:end, :) = Y(K:end, :) + b(k) * X(K-k+1:end-k+1, :);
end

n = size(X, 1);
full = conv(double(rows.valid), ones(K, 1));  % = K when all K pulses are valid
rows.valid = full(1:n) == K;
end
