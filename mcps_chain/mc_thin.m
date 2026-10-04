function [Y, xCells, yRows] = mc_thin(X, maxRows, maxCells)
%MC_THIN Thin an image for the screen: maximum of each block of rows / cells.
%
%   [Y, xCells, yRows] = mc_thin(X, maxRows, maxCells)
%
%   Y keeps at most maxRows x maxCells values, so a single detection stays
%   visible. xCells / yRows are the first cell / row of each block (for imagesc).

[Y, kr] = blockMax(X, maxRows);
[Y, kc] = blockMax(Y.', maxCells);
Y = Y.';
xCells = (0:size(Y, 2) - 1) * kc + 1;
yRows  = (0:size(Y, 1) - 1) * kr + 1;
end

function [Y, k] = blockMax(X, nMax)
% max over blocks of k rows, so that at most nMax rows remain
k = max(ceil(size(X, 1) / nMax), 1);
if k == 1
    Y = X;
    return
end
n = ceil(size(X, 1) / k) * k;
X(end + 1:n, :) = 0;
Y = reshape(max(reshape(X, k, n / k * size(X, 2)), [], 1), n / k, size(X, 2));
end
