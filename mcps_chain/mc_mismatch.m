function h = mc_mismatch(code, nTaps, spc)
%MC_MISMATCH Least-squares mismatched filter of a phase code (used by the test signal).
%
%   h = mc_mismatch(code, nTaps, spc)
%
%   code  : chip values (Nc chips)
%   nTaps : filter length in samples at fs (e.g. 210 = 3 x 35 chips x 2 samples)
%   spc   : samples per chip (fs / bandwidth = 2)
%   h     : FIR taps at fs (coefForm 'fir'). The filter is designed per chip
%           (conv(code, hc) as close as possible to one peak in its middle) and
%           every tap is held for spc samples, like the code. The peak lag is
%           (N + nTaps)/2 - 1 samples (N = code samples): the filter is centred.

c  = double(code(:));
Nc = numel(c);
Lc = round(nTaps / spc);
C  = toeplitz([c; zeros(Lc - 1, 1)], [c(1), zeros(1, Lc - 1)]);   % conv(c, hc) = C * hc
d  = zeros(Nc + Lc - 1, 1);
d(round((Nc + Lc) / 2)) = 1;
hc = C \ d;
h  = hc(min(floor((0:nTaps - 1)' / spc) + 1, Lc));                 % hold each tap spc samples
end
