function [S, idx] = rsp_nci(D, C, idx)
%RSP_NCI Non-coherent integration of every Doppler bin over n consecutive frames.
%
%   [S, idx] = rsp_nci(D, C, idx)
%
%   D   : frames x range x doppler (complex FFT output)
%   C   : P.nci (see rsp_default_params)
%   idx : pulse index of each frame (default 1:size(D,1))
%   S   : frames x range x doppler (real), for every range cell and bin
%           S(f) = sum_{i=0}^{n-1} law(D(f-i))       (or the mean if C.average)
%         law = |.|^2 ('square'), |.| ('linear') or 20log10|.| ('log')
%   idx : pulse index of each row of S (the newest frame in the buffer)

nF = size(D, 1);
if nargin < 3 || isempty(idx)
    idx = (1:nF)';
end
idx = idx(:);
n   = C.nFrames;
if n > nF
    error('rsp_nci:frames', 'Buffer of %d frames, only %d available.', n, nF);
end

switch lower(C.law)
    case 'square', S = real(D).^2 + imag(D).^2;
    case 'linear', S = abs(D);
    case 'log',    S = 20*log10(abs(D) + eps(class(D)));
    otherwise, error('rsp_nci:law', 'Unknown P.nci.law "%s".', C.law);
end

if n > 1
    S = filter(ones(n, 1, class(S)), 1, S, [], 1);       % running sum over frames
end
if C.average
    S = S / n;
end

switch lower(C.output)
    case 'valid'
        S   = S(n:end, :, :);
        idx = idx(n:end);
    case 'same'
    otherwise
        error('rsp_nci:output', 'Unknown P.nci.output "%s".', C.output);
end
end
