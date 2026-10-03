function [D, idx, w] = rsp_doppler_fft(X, F, idx)
%RSP_DOPPLER_FFT Slow-time FFT of every range cell (sliding or block CPI).
%
%   [D, idx, w] = rsp_doppler_fft(X, F, idx)
%
%   X   : pulses x range (complex)
%   F   : P.fft (see rsp_default_params)
%   idx : pulse index of each row of X (default 1:size(X,1))
%   D   : frames x range x nfft (complex), D(f,r,b) = Doppler bin b of the
%         FFT that ends at pulse idx(f)
%   idx : pulse index of the LAST pulse of each FFT window
%   w   : slow-time window used
%
%   Frame f uses pulses (f-1)*hop + (1:nPulses). Frames are processed in
%   blocks of F.chunk to keep the temporary memory small.

[M, R] = size(X);
if nargin < 3 || isempty(idx)
    idx = (1:M)';
end
idx  = idx(:);
N    = F.nPulses;
nfft = max(F.nfft, N);
hop  = F.hop;
cls  = class(X);

if M < N
    error('rsp_doppler_fft:pulses', 'Only %d pulses, the FFT needs %d.', M, N);
end

w = rsp_window(F.window, N, F.windowParam);
switch lower(F.norm)
    case 'none'
    case 'noise', w = w / norm(w);      % noise power gain = 1
    case 'peak',  w = w / sum(w);       % unit gain for a tone on a bin centre
    otherwise, error('rsp_doppler_fft:norm', 'Unknown P.fft.norm "%s".', F.norm);
end
w = cast(w, cls);

starts = 1:hop:(M - N + 1);
nF     = numel(starts);
D      = complex(zeros(nF, R, nfft, cls));
chunk  = max(1, F.chunk);

for c = 1:chunk:nF
    f    = c:min(c + chunk - 1, nF);
    rows = bsxfun(@plus, starts(f), (0:N-1)');          % N x numel(f)
    blk  = reshape(X(rows(:), :), N, numel(f), R);      % N x frames x range
    blk  = blk .* w;                                    % slow-time window
    S    = fft(blk, nfft, 1);                           % nfft x frames x range
    if F.shift
        S = fftshift(S, 1);
    end
    D(f, :, :) = permute(S, [2 3 1]);
end

idx = idx(starts + N - 1);
end
