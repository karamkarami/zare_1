function P = mc_params()
%MC_PARAMS Every parameter of the chain. Edit the values in this file.
%
%   video -> decoder -> 3-pulse canceler -> Doppler FFT -> NCI buffer -> SO-CFAR -> max
%
%   Data layout (same as the MCPS log):
%     video, decoder, canceler : pulses x range cells            (complex)
%     integral, cfar           : rows x range cells x Doppler bins (real)
%   Units: MHz and microseconds (MHz * us = number of cycles / samples).

%% Sampling ------------------------------------------------------------------------------
P.fsMHz = 2;                 % complex sampling rate [MHz] -> range cell = c / (2 fs) = 75 m
                             % one PRF is about 5469 samples (= range cells)

%% Pulses --------------------------------------------------------------------------------
% Pulse 1 = short pulse, pulse 2 = long pulse. Both are phase codes of 1 MHz bandwidth
% (1 us chips = 2 samples per chip at 2 MHz). Each one has its own mismatched filter.
%
%   coef       decoder (mismatched filter) coefficients at fs, one value per sample.
%              Length 210 (short) and 1440 (long) = 3 x pulse length in samples.
%              [] = placeholder filter (the chain runs, but this is NOT your decoder).
%   coefFile   instead of coef: a .mat file (first numeric variable) or a text file
%              (one value per line, or two columns "real imag").
%   code       OPTIONAL: the transmitted chips (35 / 240 values, +-1 or complex).
%              Only used to check the decoder (peak position, PSL, loss) and by the
%              synthetic test. Not needed to process the data.
%   freqMHz    centre frequency of the pulse inside the video [MHz]. 0 = both pulses
%              on the same frequency. To separate them later, e.g. -0.5 and +0.5, and
%              set P.decoder.bandFilter = true.
%   txDelayUs  transmit start of this pulse relative to range cell 1 [us]. The echo of
%              range cell r starts at sample r + txDelay. If the short pulse is sent
%              right after the long one, its delay is the long pulse length (240 us).
%              main_chain prints the delay that best fits the log decoder.
%   alignLag   [] = automatic (peak of the filter), or the peak lag in samples.
%   gain       linear gain of this pulse's output (to level the two parts).

P.pulse(1).name       = 'short';
P.pulse(1).widthUs    = 35;          % 35 us -> 70 samples
P.pulse(1).bwMHz      = 1;
P.pulse(1).coef       = [];          % <<< put the 210 coefficients of the short decoder here
P.pulse(1).coefFile   = '';          %     ... or the file that holds them
P.pulse(1).coefLength = 210;         % expected number of coefficients (checked)
P.pulse(1).code       = [];          % optional: the 35 chips of the short pulse
P.pulse(1).freqMHz    = 0;
P.pulse(1).txDelayUs  = 240;         % sent after the long pulse (check with main_chain)
P.pulse(1).alignLag   = [];
P.pulse(1).gain       = 1;

P.pulse(2).name       = 'long';
P.pulse(2).widthUs    = 240;         % 240 us -> 480 samples
P.pulse(2).bwMHz      = 1;
P.pulse(2).coef       = [];          % <<< put the 1440 coefficients of the long decoder here
P.pulse(2).coefFile   = '';
P.pulse(2).coefLength = 1440;
P.pulse(2).code       = [];          % optional: the 240 chips of the long pulse
P.pulse(2).freqMHz    = 0;
P.pulse(2).txDelayUs  = 0;           % sent first
P.pulse(2).alignLag   = [];
P.pulse(2).gain       = 1;

%% Decoder --------------------------------------------------------------------------------
P.decoder.switchCell  = 560;         % cells 1..559 from the short pulse, 560..end from the long one
P.decoder.coefForm    = 'fir';       % 'fir'       : coefficients are the impulse response h,
                                     %               y = conv(x, h)
                                     % 'correlate' : coefficients are a reference c,
                                     %               y(n) = sum_m x(n+m) conj(c(m))
P.decoder.normalize   = 'none';      % 'none' | 'noise' (divide by norm(h): same noise
                                     % floor for both pulses)
P.decoder.bandFilter  = false;       % low-pass filter after moving each pulse to 0 Hz
                                     % (turn on when the two frequencies are different)
P.decoder.bandFactor  = 1;           % low-pass cut-off = bandFactor * bwMHz / 2
P.decoder.bandTaps    = 41;          % low-pass length (odd)
P.decoder.blockPulses = 256;         % pulses per FFT block (memory only)

%% 3-pulse canceler ------------------------------------------------------------------------
P.canceler.enable = true;
P.canceler.coef   = [1 -2 1];        % y(n) = x(n) - 2 x(n-1) + x(n-2)   (2-pulse: [1 -1])

%% Doppler FFT ------------------------------------------------------------------------------
P.doppler.nPulses  = 16;             % pulses in one FFT
P.doppler.nfft     = 16;             % FFT size = number of Doppler bins (>= nPulses)
P.doppler.hop      = 1;              % pulses between two FFTs: 1 = one FFT per pulse
                                     % (sliding), 16 = consecutive blocks
P.doppler.window   = 'rect';         % 'rect' | 'hann' | 'hamming' | 'blackman' | vector
P.doppler.fftshift = false;          % true = zero Doppler in the middle bin
P.doppler.detector = 'abs';          % 'abs' |X| or 'power' |X|^2 (input of the integration)
P.doppler.alignToLog    = true;      % compute the FFTs that end on the same pulses as the
                                     % rows of the log "integral" lane (when it is loaded)
P.doppler.logPulseShift = 0;         % pulses to add to the log row label (0 = the log
                                     % labels a row with its newest pulse)

%% Non-coherent integration ------------------------------------------------------------------
P.nci.length = 4;                    % buffer length: consecutive FFT outputs summed per bin
                                     % (1 = no integration)
P.nci.mode   = 'sum';                % 'sum' | 'mean'

%% CFAR (along range, for every Doppler bin) ------------------------------------------------
P.cfar.type        = 'SO';           % 'SO' smallest-of | 'CA' cell-averaging | 'GO' greatest-of
P.cfar.nRef        = 16;             % reference cells on EACH side
P.cfar.nGuard      = 4;              % guard cells on EACH side
P.cfar.thresholdDb = 15;             % threshold above the noise estimate [dB, power]
P.cfar.minNoise    = 0;              % lower limit of the noise estimate (same units as
                                     % the integration output; 0 = none)
P.cfar.edge        = 'partial';      % 'partial' : near the ends use the cells that exist
                                     % 'skip'    : no detection without both full windows
P.cfar.cells       = [];             % [first last] cells tested, [] = all (e.g. [70 Inf]
                                     % leaves out the eclipsed cells of the short pulse)
P.cfar.output      = 'value';        % value kept on a detection (all other cells 0):
                                     % 'value' (integrated value) | 'ratio' (value / noise
                                     % estimate) | 'binary' (1)

%% Output: max over Doppler bins --------------------------------------------------------------
P.output.bins = [];                  % bins entering the max: [] = all, e.g. 2:16 = no zero Doppler

%% Comparison with the log --------------------------------------------------------------------
P.compare.stageWise = true;          % also feed each block with the previous LOG lane, so a
                                     % difference points at one block
P.compare.maxShift  = 1000;          % range shift searched [cells]
P.compare.maxRows   = 600;           % rows used for the numbers (memory)
P.compare.plot      = true;          % figures: mine / log / difference / range profile
end
