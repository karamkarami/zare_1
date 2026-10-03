function out = rsp_cfar_cells(S, P, cells, nInt, nRef)
%RSP_CFAR_CELLS CFAR on a range interval only; the rest is "not tested".
%
%   out = rsp_cfar_cells(S, P, cells, nInt, nRef)
%
%   S     : frames x range x doppler
%   cells : [first last] range cells to process. Cells outside (partly
%           received echoes) are neither tested nor used as reference, so
%           their lower noise cannot pull a smallest-of estimate down.
%   Output as rsp_cfar, full size, with range indices of the full array.

[nF, R, nB] = size(S);
c   = cells(1):cells(2);
sub = rsp_cfar(S(:, c, :), P.cfar, P.nci.law, nInt, nRef);

out.det       = false(nF, R, nB);
out.threshold = inf(nF, R, nB, 'single');
out.map       = zeros(nF, R, nB, 'single');
out.det(:, c, :)       = sub.det;
out.threshold(:, c, :) = sub.threshold;
out.map(:, c, :)       = sub.map;
out.factor = sub.factor;
out.list   = sub.list;
out.list.range = out.list.range + c(1) - 1;
end
