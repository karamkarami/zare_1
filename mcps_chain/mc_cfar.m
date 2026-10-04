function [D, T] = mc_cfar(Y, P)
%MC_CFAR CFAR along range, separately for every row and every Doppler bin.
%
%   [D, T] = mc_cfar(Y, P)
%
%   Y : rows x range cells x bins (integration output)
%   D : same size: the detected cells keep a value (P.cfar.output), all others are 0
%   T : threshold, same size (only computed when asked for)
%
%   Cell under test r, with G guard and R reference cells on each side:
%       lead  = mean of Y(r-G-R .. r-G-1)          cells before
%       lag   = mean of Y(r+G+1 .. r+G+R)          cells after
%       noise = min(lead, lag)                     SO (CA: mean of both, GO: max)
%       T     = alpha * max(noise, minNoise)       detection when Y(r) > T
%   alpha = 10^(thresholdDb/20) for an amplitude |X| and 10^(thresholdDb/10)
%   for a power |X|^2, so thresholdDb is always a power ratio.
%   Near the ends a window has fewer cells ('partial': mean of the cells that
%   exist, one side only when the other is empty; 'skip': no detection).
%   Only the cells P.cfar.cells = [first last] are tested ([] = all).

C = P.cfar;
[nRows, nCells, nBins] = size(Y);
if strcmpi(P.doppler.detector, 'power')
    alpha = 10^(C.thresholdDb / 10);
else
    alpha = 10^(C.thresholdDb / 20);
end
minNoise = 0;
if isfield(C, 'minNoise') && ~isempty(C.minNoise)
    minNoise = C.minNoise;
end

% window limits of every cell (the same for all rows and bins)
r  = 1:nCells;
G  = C.nGuard;
R  = C.nRef;
a1 = max(r - G - R, 1);       b1 = r - G - 1;                 % lead window a1..b1
a2 = r + G + 1;               b2 = min(r + G + R, nCells);    % lag window  a2..b2
n1 = max(b1 - a1 + 1, 0);     n2 = max(b2 - a2 + 1, 0);       % cells in each window
b1 = max(b1, a1 - 1);                                         % empty window -> sum 0
a2 = min(a2, nCells + 1);     b2 = max(b2, a2 - 1);
test = true(1, nCells);                                       % cells that are tested
if strcmpi(C.edge, 'skip')
    test = n1 == R & n2 == R;
end
if isfield(C, 'cells') && ~isempty(C.cells)
    test = test & r >= C.cells(1) & r <= C.cells(2);
end

D = zeros(nRows, nCells, nBins, 'single');
if nargout > 1
    T = D;
end
for k = 1:nBins
    A  = double(Y(:, :, k));
    CS = [zeros(nRows, 1), cumsum(A, 2)];          % CS(:, j+1) = sum of cells 1..j
    s1 = CS(:, b1 + 1) - CS(:, a1);                % sum of the lead window
    s2 = CS(:, b2 + 1) - CS(:, a2);                % sum of the lag window
    m1 = s1 ./ n1;                                 % mean (NaN when the window is empty;
    m2 = s2 ./ n2;                                 % min / max ignore NaN)
    switch upper(C.type)
        case 'SO', noise = min(m1, m2);
        case 'GO', noise = max(m1, m2);
        case 'CA', noise = (s1 + s2) ./ (n1 + n2);
        otherwise, error('mc_cfar:type', 'Unknown P.cfar.type "%s".', C.type);
    end
    noise(noise < minNoise) = minNoise;            % NaN stays NaN -> no detection
    thr   = alpha * noise;
    det   = A > thr & test;

    out = zeros(nRows, nCells);
    switch lower(C.output)
        case 'value',  out(det) = A(det);
        case 'ratio',  out(det) = A(det) ./ noise(det);
        case 'binary', out(det) = 1;
        otherwise, error('mc_cfar:output', 'Unknown P.cfar.output "%s".', C.output);
    end
    D(:, :, k) = out;
    if nargout > 1
        T(:, :, k) = thr;
    end
end
end
