function h = mc_lowpass(cutoffMHz, fsMHz, nTaps)
%MC_LOWPASS Linear-phase low-pass FIR (Hamming-windowed sinc), unit gain at 0 Hz.
%
%   h = mc_lowpass(cutoffMHz, fsMHz, nTaps)     cut-off = -6 dB edge [MHz]
%
%   The delay of the filter is (nTaps - 1) / 2 samples (nTaps is made odd).

nTaps = nTaps + mod(nTaps + 1, 2);                   % odd -> integer delay
n  = (0:nTaps-1)' - (nTaps - 1) / 2;
fc = cutoffMHz / fsMHz;                              % cycles per sample
x  = 2 * fc * n;
h  = ones(nTaps, 1);
h(x ~= 0) = sin(pi * x(x ~= 0)) ./ (pi * x(x ~= 0));
h  = h .* mc_window('hamming', nTaps);
h  = h / sum(h);
end
