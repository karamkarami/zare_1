function [S, rows] = mc_doppler(X, P, rows, endPulses)
%MC_DOPPLER Doppler FFT along slow time, for every range cell.
%
%   [S, rows] = mc_doppler(X, P, rows)              one FFT every P.doppler.hop pulses
%   [S, rows] = mc_doppler(X, P, rows, endPulses)   FFTs that end on these rows of X
%
%   X : pulses x range cells (canceler output)
%   S : FFTs x range cells x nfft bins, detected (|F| or |F|^2), single
%
%   FFT m uses the nPulses pulses that end on pulse e = endPulses(m), oldest first:
%       F(k, r) = sum_{n=1..nPulses} w(n) X(e - nPulses + n, r) exp(-j 2 pi (k-1)(n-1) / nfft)
%   Bin k holds the Doppler (k-1)/nfft cycles per pulse (0 Hz in bin 1).

D = P.doppler;
[nPulses, nCells] = size(X);
N = D.nPulses;
if D.nfft < N
    error('mc_doppler:nfft', 'P.doppler.nfft (%d) must be >= nPulses (%d).', D.nfft, N);
end
if nargin < 4 || isempty(endPulses)
    endPulses = (N:D.hop:nPulses)';
end
endPulses = endPulses(:);
endPulses = endPulses(endPulses >= N & endPulses <= nPulses);

w = cast(mc_window(D.window, N), class(X));
S = zeros(numel(endPulses), nCells, D.nfft, 'single');
for m = 1:numel(endPulses)
    block = X(endPulses(m)-N+1 : endPulses(m), :);    % N pulses x range cells
    F = fft(block .* w, D.nfft, 1);                   % nfft bins x range cells
    if D.fftshift
        F = fftshift(F, 1);
    end
    if strcmpi(D.detector, 'power')
        F = abs(F).^2;
    else
        F = abs(F);
    end
    S(m, :, :) = reshape(F.', 1, nCells, D.nfft);
end

full = conv(double(rows.valid), ones(N, 1));          % = N when all N pulses are valid
full = full(1:nPulses) == N;
rows.valid = full(endPulses);
rows.pulse = rows.pulse(endPulses);
rows.seq   = rows.seq(endPulses);
end
