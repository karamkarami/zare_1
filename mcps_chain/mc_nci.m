function [Y, rows] = mc_nci(S, P, rows)
%MC_NCI Non-coherent integration with a buffer on every Doppler bin.
%
%   [Y, rows] = mc_nci(S, P, rows)
%
%   S : FFTs x range cells x bins (detected FFT output)
%   Y : same size. Every (range cell, bin) has a FIFO buffer that holds its
%       last L = P.nci.length FFT outputs; the output is their sum (or mean):
%           Y(m, r, k) = S(m, r, k) + S(m-1, r, k) + ... + S(m-L+1, r, k)
%       The first L-1 rows are partial sums and are marked rows.valid = false.

L = P.nci.length;
[nRows, nCells, nBins] = size(S);
Y      = zeros(nRows, nCells, nBins, 'single');
buffer = zeros(L, nCells, nBins, 'single');   % last L outputs of every cell and bin
slot   = 1;                                   % slot of the oldest output
for m = 1:nRows
    buffer(slot, :, :) = S(m, :, :);          % the newest output replaces the oldest
    slot = mod(slot, L) + 1;
    Y(m, :, :) = sum(buffer, 1);
end
if strcmpi(P.nci.mode, 'mean')
    Y = Y / L;
end

full = conv(double(rows.valid), ones(L, 1));  % = L when all L outputs are valid
rows.valid = full(1:nRows) == L;
end
