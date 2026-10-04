function [Y, each, info] = mc_decoder(X, P)
%MC_DECODER Pulse compression of the two phase-coded pulses (mismatched filters).
%
%   [Y, each, info] = mc_decoder(X, P)
%
%   X    : video, pulses x range samples (complex)
%   Y    : decoder output, pulses x range cells. Cells 1 .. switchCell-1 come from
%          pulse 1 (short), cells switchCell .. end from pulse 2 (long).
%   each : {output of pulse 1, output of pulse 2}, each over the whole range
%          (computed only when asked for; used to check the alignment)
%   info : one entry per pulse, see mc_filter
%
%   For every pulse k:
%     1. x_k = X .* exp(-j 2 pi f_k n / fs)   move the pulse to 0 Hz (nothing when f_k = 0)
%     2. z_k = x_k (*) h_k                     FIR filter along range (FFT convolution);
%                                              h_k = [low-pass (*)] mismatched filter
%     3. Y_k(r) = z_k(r + shift_k)             shift_k = peak lag of the filter + transmit
%                                              delay, so the peak of a target in range
%                                              cell r lands in cell r

[nPulses, nCells] = size(X);
nType    = numel(P.pulse);
wantEach = nargout > 1;

for k = 1:nType
    info(k) = mc_filter(P.pulse(k), P);  %#ok<AGROW>
end

% which pulse feeds which output cell
source = ones(1, nCells);
if nType > 1
    source(min(P.decoder.switchCell, nCells + 1):end) = 2;
end

Lmax = max(arrayfun(@(s) numel(s.h), info));
nfft = 2^nextpow2(nCells + Lmax - 1);           % no circular wrap-around
n    = 0:nCells-1;
blk  = P.decoder.blockPulses;

Y    = complex(zeros(nPulses, nCells, class(X)));
each = cell(1, nType);
for k = 1:nType
    cellsK = find(source == k);
    if isempty(cellsK) && ~wantEach
        continue
    end
    if wantEach
        each{k} = complex(zeros(nPulses, nCells, class(X)));
    end
    h   = info(k).h;
    H   = cast(fft(h, nfft).', class(X));                       % 1 x nfft
    f   = P.pulse(k).freqMHz / P.fsMHz;                         % cycles per sample
    mix = cast(exp(-1j*2*pi*f*n), class(X));                    % 1 x nCells
    idx = (1:nCells) + info(k).shift;                           % sample of z read by each cell
    ok  = idx >= 1 & idx <= nCells + numel(h) - 1;

    for first = 1:blk:nPulses
        rows = first:min(first + blk - 1, nPulses);
        x = X(rows, :);
        if f ~= 0
            x = x .* mix;                                       % 1. pulse to 0 Hz
        end
        z  = ifft(fft(x, nfft, 2) .* H, [], 2);                 % 2. linear convolution
        yk = complex(zeros(numel(rows), nCells, class(X)));
        yk(:, ok) = z(:, idx(ok));                              % 3. range alignment
        Y(rows, cellsK) = yk(:, cellsK);
        if wantEach
            each{k}(rows, :) = yk;
        end
    end
end
end
