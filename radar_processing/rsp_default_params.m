function P = rsp_default_params()
%RSP_DEFAULT_PARAMS Default parameter set of the radar signal processing chain.
%
%   P = rsp_default_params() returns a struct holding every tunable value:
%
%       video -> decoder (matched / mismatched filter) -> 3-pulse canceler
%             -> Doppler FFT -> non-coherent integration -> SO-CFAR
%             -> max over Doppler bins -> plots / PPI
%
%   Edit the returned struct (or this file) to match your radar.
%
%   Data layout everywhere in the chain:
%       2-D lanes : pulses x range            (rows = slow time, cols = fast time)
%       3-D lanes : pulses x range x doppler  (same layout as the MCPS log)
%   Units: MHz and microseconds inside the signal chain (MHz * us = cycles),
%          metres, degrees, m/s and Hz for the radar geometry.

%% ------------------------------------------------------------------ Radar
P.radar.prfHz        = 1000;        % pulse repetition frequency [Hz]
P.radar.rpm          = 6;           % antenna rotation [rev/min] (36 deg/s)
P.radar.fcMHz        = 600;         % carrier frequency [MHz] (UHF, lambda = 0.5 m):
                                    % unambiguous velocity +-lambda*PRF/4 = +-125 m/s
P.radar.azStartDeg   = 0;           % antenna azimuth of pulse 1 [deg], 0 = north
P.radar.rangeOffsetM = 0;           % range of cell 1 [m]

%% ------------------------------------------------------------------ Antenna (azimuth only)
P.antenna.beamwidthDeg = 4;         % 3 dB beamwidth [deg]
P.antenna.sidelobeDb   = -25;       % first sidelobe, one way [dB]
P.antenna.backlobeDb   = -40;       % back lobe, one way [dB]
P.antenna.nullDepthDb  = -20;       % filling of the nulls between sidelobes [dB]

%% ------------------------------------------------------------------ Sampling
P.fs = 6;                           % complex baseband sampling frequency [MHz]
                                    % (6 samples per 1 us chip, 25 m range cells)

%% ------------------------------------------------------------------ Pulses
% One entry per transmitted pulse. The echo of a scatterer in range cell r
% (0-based) starts at sample r + delay of that pulse.
%
%   type         'code' (phase code) | 'lfm'
%   bwMHz        bandwidth [MHz]: chip rate of a code (1 MHz -> 1 us chips),
%                swept bandwidth of an LFM
%   code         chip values of the code, any length, binary or polyphase
%                (put your own coefficients here; rsp_code has common codes)
%   codeRate     'chip' (one value per chip) | 'sample' (already sampled at fs)
%   decoder      your decoder coefficients; [] = matched filter
%   decoderForm  'fir'       : FIR taps applied by convolution
%                              (matched filter = conj(fliplr(code)))
%                'reference' : correlation reference (matched filter = code)
%   decoderRate  'chip' | 'sample'
%   decoderLag   [] = align on the main peak automatically, or lag in samples
%   widthUs      LFM length [us] (a code's length is numel(code)/bwMHz)
%   slope        LFM: +1 up-chirp, -1 down-chirp
%   fcMHz        pulse centre frequency inside the baseband [MHz]
%   delayUs      transmit start of this pulse relative to range cell 0 [us]
%   window       weighting of the decoder ('none', 'hamming', 'taylor' ...)
%   windowParam  kaiser -> beta, taylor -> [nbar sllDb]
%   gainDb       extra gain of this pulse's decoder output
%   samples      user transmit samples at fs (overrides type / code / LFM)

% Pulse 1: short pulse, Barker 13 (13 us) with a 52-tap mismatched decoder
% (PSL -48 dB instead of -22 dB, 0.2 dB loss); covers the blind zone of pulse 2
P.pulse(1).type        = 'code';
P.pulse(1).bwMHz       = 1;
P.pulse(1).code        = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
P.pulse(1).codeRate    = 'chip';
P.pulse(1).decoder     = rsp_code_mmf(P.pulse(1).code, 52);   % [] = matched filter
P.pulse(1).decoderForm = 'fir';
P.pulse(1).decoderRate = 'chip';
P.pulse(1).decoderLag  = [];
P.pulse(1).widthUs     = 13;
P.pulse(1).slope       = 1;
P.pulse(1).fcMHz       = -1.5;      % frequency diversity: short and long pulses
                                    % apart, so the long echo does not leak into
                                    % the short decoder (set 0 if your radar uses one frequency)
P.pulse(1).delayUs     = 100;       % sent right after pulse 2
P.pulse(1).window      = 'none';
P.pulse(1).windowParam = [];
P.pulse(1).gainDb      = 0;
P.pulse(1).samples     = [];

% Pulse 2: long pulse, LFM 100 us, 1 MHz, Taylor weighted (PSL -41 dB)
P.pulse(2)             = P.pulse(1);
P.pulse(2).type        = 'lfm';
P.pulse(2).code        = [];
P.pulse(2).decoder     = [];
P.pulse(2).widthUs     = 100;
P.pulse(2).bwMHz       = 1;
P.pulse(2).slope       = 1;
P.pulse(2).delayUs     = 0;         % sent first
P.pulse(2).fcMHz       = 1.5;
P.pulse(2).window      = 'taylor';
P.pulse(2).windowParam = [6 -45];

% Phase code instead of the LFM (e.g. a 63-chip m-sequence):
%   P.pulse(2).type = 'code';  P.pulse(2).code = rsp_code('mls', 63);
%   P.pulse(2).window = 'none';  P.pulse(1).delayUs = 63;

%% ------------------------------------------------------------------ Decoder
P.mf.enable      = true;
P.mf.combine     = 'stitch';    % 'stitch' | 'pulse1' | 'pulse2'
P.mf.switchCell  = [];          % cells 1..switchCell from the short pulse
                                % [] = automatic = blind zone of the long pulse
P.mf.norm        = 'noise';     % 'noise' : unit noise gain (same floor for both pulses)
                                % 'peak'  : unit gain for a matched echo
                                % 'none'  : raw
P.mf.rangeOffset = 0;           % extra system delay to remove [samples]
P.mf.channelFilter   = true;    % band-pass each pulse's frequency channel before
                                % its decoder (dual-frequency receiver); keeps the
                                % strong long echo out of the short decoder.
                                % false = plain decoder (e.g. to match a log)
P.mf.channelBwFactor = 3;       % passband = factor * pulse bandwidth (3: keeps the
                                % Barker mismatched decoder at PSL -46 dB)
P.mf.channelRejectDb = 60;      % stopband rejection [dB]
P.mf.keepEach    = false;       % also return each pulse's own output

%% ------------------------------------------------------------------ Canceler (MTI)
P.canceler.enable    = true;
P.canceler.order     = 3;       % number of pulses: 2, 3 (default), 4 ...
P.canceler.weights   = [];      % custom slow-time taps, e.g. [1 -2 1]; [] = binomial
P.canceler.normalize = false;   % true = unit noise gain
P.canceler.output    = 'same';  % 'same' (first rows 0) | 'valid' (drop them)

%% ------------------------------------------------------------------ Doppler FFT
P.fft.nPulses     = 16;         % pulses per FFT
P.fft.nfft        = 16;         % FFT size (>= nPulses)
P.fft.hop         = 1;          % pulses between FFTs (1 = one FFT per pulse)
P.fft.window      = 'hamming';
P.fft.windowParam = [];
P.fft.norm        = 'none';     % 'none' | 'noise' | 'peak'
P.fft.shift       = false;      % true = zero Doppler in the middle bin
P.fft.chunk       = 64;         % frames per block (memory control)

%% ------------------------------------------------------------------ Non-coherent integration
P.nci.nFrames = 32;             % n : consecutive FFT outputs per bin. With hop = 1
                                % the integration spans nPulses + n - 1 = 47 pulses
                                % (1.7 deg at 6 rpm), inside the 111-pulse beam dwell;
                                % more looks -> lower CFAR threshold
P.nci.law     = 'square';       % 'square' |x|^2 | 'linear' |x| | 'log' dB
P.nci.average = true;           % mean (true) or sum (false)
P.nci.output  = 'valid';        % 'valid' | 'same'

%% ------------------------------------------------------------------ CFAR
P.cfar.type          = 'SO';    % 'SO' | 'CA' | 'GO'
P.cfar.nRef          = [];      % reference cells on EACH side
                                % [] = 16 resolution cells (rsp_cfar_window)
P.cfar.nGuard        = [];      % guard cells on EACH side
                                % [] = widest compressed main lobe + 2
P.cfar.thresholdMode = 'pfa';   % 'pfa' | 'factor'
P.cfar.pfa           = 1e-6;
P.cfar.factorDb      = 13;      % used when thresholdMode = 'factor'
P.cfar.nIntEffective = [];      % [] = from canceler, window, hop, nFrames (rsp_cfar_looks)
P.cfar.nRefEffective = [];      % [] = from the decoder's noise correlation
P.cfar.edge          = 'oneSided';  % 'oneSided' | 'partial' | 'none'
P.cfar.bins          = [];      % Doppler bins tested ([] = all)
P.cfar.rangeCells    = [];      % [first last] cells tested ([] = all)
P.cfar.mapValue      = 'value'; % CFAR output on detections: 'value' (integrated
                                % value, like the log lane) | 'snr' (value / noise
                                % estimate: bins compare fairly in the max stage)
P.cfar.validOnly     = true;    % run only on cells whose echo is received in full
                                % (outside: partly eclipsed or cut by the record end)

%% ------------------------------------------------------------------ Output (after CFAR)
P.output.bins        = [];      % bins taken into the max after CFAR ([] = all),
                                % e.g. 5 to see bin 5 only, or [2:16] without zero Doppler
P.output.keepDoppler = true;    % keep the complex FFT output in rsp_chain

%% ------------------------------------------------------------------ Plot extraction
P.plots.gapPulses = 3;          % detections closer than this (pulses) ...
P.plots.gapCells  = 6;          % ... and this (range cells) belong to one plot
P.plots.minHits     = 10;       % fewer detections than this = discarded
P.plots.minWidthDeg = 1;        % narrower in azimuth than this = discarded
                                % (a target lasts about one beamwidth)

%% ------------------------------------------------------------------ PPI display
P.ppi.source      = 'max';      % 'max' (output after CFAR) | 'video' | 'decoder' | 'canceler'
P.ppi.maxScale    = 'margin';   % 'max' painted as: 'margin' (dB over the CFAR
                                % threshold) | 'value' (integrated value, dB)
P.ppi.maxRangeKm  = [];         % [] = whole record
P.ppi.pixels      = 800;        % image size
P.ppi.ringKm      = 20;         % range ring spacing
P.ppi.climDb      = [];         % colour limits [dB], [] = automatic
P.ppi.fadeDb      = 12;         % afterglow: fading over one revolution [dB]
P.ppi.showPlots   = true;       % mark extracted plots
P.ppi.trailScans  = 4;          % plots of the last revolutions kept as a trail
P.ppi.colormap    = 'phosphor'; % 'phosphor' or any MATLAB colormap name
P.ppi.blockPulses = 250;        % pulses processed per block
P.ppi.stepDeg     = 0.5;        % sweep step per screen update [deg]
P.ppi.speed       = 1;          % sweep speed: 1 = real antenna speed, Inf = as fast
                                % as possible (processing may be slower)
end
