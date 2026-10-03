function P = rsp_default_params()
%RSP_DEFAULT_PARAMS Default parameter set of the radar signal processing chain.
%
%   P = rsp_default_params() returns a struct holding every tunable value of
%   the chain:
%
%       video -> matched filter -> 3-pulse canceler -> Doppler FFT
%             -> non-coherent integration -> SO-CFAR
%
%   Edit the returned struct (or this file) to match your radar. All values
%   below are only "reasonable" defaults used by the self-test.
%
%   Data convention used everywhere in the chain:
%       2-D lanes : pulses x range            (rows = slow time, cols = fast time)
%       3-D lanes : pulses x range x doppler  (same layout as the MCPS log)
%
%   Units: frequencies in MHz, times in microseconds (MHz * us = cycles).

%% ------------------------------------------------------------------ Sampling
P.fs = 4;                       % complex baseband sampling frequency [MHz]

%% ------------------------------------------------------------------ Pulses
% One entry per transmitted pulse. The echo of a scatterer located at range
% cell r appears in the video at samples  r + delay : r + delay + width - 1.
%
%   widthUs     pulse length [us]                    (p1, p2)
%   bwMHz       LFM swept bandwidth [MHz]            (0 = unmodulated pulse)
%   fcMHz       pulse centre frequency inside the baseband [MHz]
%   slope       +1 up-chirp, -1 down-chirp
%   delayUs     transmit start of this pulse, relative to range cell 0 [us]
%   code        optional phase code (chip values, e.g. Barker [1 1 1 -1 1]);
%               chips are spread uniformly over widthUs. [] = not used
%   samples     optional user replica (complex column vector at fs).
%               When not empty it overrides every field above except delayUs
%   window      matched-filter amplitude weighting:
%               'none' | 'hamming' | 'hann' | 'blackman' | 'blackmanharris'
%               | 'kaiser' | 'taylor' | numeric vector | function handle @(N)
%   windowParam window parameter: kaiser -> beta, taylor -> [nbar sllDb]
%   gainDb      extra gain applied to this pulse's matched-filter output

% Pulse 1: short pulse (covers the blind zone of the long pulse)
P.pulse(1).widthUs     = 2;         % p1
P.pulse(1).bwMHz       = 0;
P.pulse(1).fcMHz       = -1.2;
P.pulse(1).slope       = +1;
P.pulse(1).delayUs     = 20;        % transmitted right after the long pulse
P.pulse(1).code        = [];
P.pulse(1).samples     = [];
P.pulse(1).window      = 'none';
P.pulse(1).windowParam = [];
P.pulse(1).gainDb      = 0;

% Pulse 2: long LFM pulse
P.pulse(2).widthUs     = 20;        % p2
P.pulse(2).bwMHz       = 2;
P.pulse(2).fcMHz       = 0.8;
P.pulse(2).slope       = +1;
P.pulse(2).delayUs     = 0;         % transmitted first
P.pulse(2).code        = [];
P.pulse(2).samples     = [];
P.pulse(2).window      = 'taylor';
P.pulse(2).windowParam = [4 -35];   % [nbar sidelobe level dB]
P.pulse(2).gainDb      = 0;

%% ------------------------------------------------------------------ Matched filter
P.mf.enable      = true;
P.mf.combine     = 'stitch';    % 'stitch' | 'pulse1' | 'pulse2' (any 'pulseK')
P.mf.switchCell  = [];          % cells 1..switchCell come from the short pulse,
                                % the rest from the long pulse.
                                % [] = automatic = blind zone of the long pulse
P.mf.norm        = 'noise';     % 'noise' : unit noise gain (noise floor kept)
                                % 'peak'  : unit amplitude gain for a matched echo
                                % 'none'  : raw correlation
P.mf.rangeOffset = 0;           % extra system delay to remove [samples]
P.mf.keepEach    = true;        % also return each pulse's own MF output

%% ------------------------------------------------------------------ Canceler (MTI)
P.canceler.enable    = true;
P.canceler.order     = 3;       % number of pulses: 2, 3 (default), 4, ...
P.canceler.weights   = [];      % custom slow-time taps (overrides order),
                                % e.g. [1 -2 1]. [] = binomial weights
P.canceler.normalize = false;   % true = divide taps by their norm (unit noise gain)
P.canceler.output    = 'same';  % 'same'  : keep all pulses, transient rows = 0
                                % 'valid' : drop the (order-1) transient rows

%% ------------------------------------------------------------------ Doppler FFT
P.fft.nPulses     = 16;         % pulses per FFT (CPI length)
P.fft.nfft        = 16;         % FFT size (>= nPulses, zero padded)
P.fft.hop         = 1;          % pulses between consecutive FFTs
                                % 1 = sliding (one output per pulse)
                                % nPulses = non-overlapping CPIs
P.fft.window      = 'hamming';  % same choices as the pulse window
P.fft.windowParam = [];
P.fft.norm        = 'none';     % 'none' | 'noise' | 'peak'
P.fft.shift       = false;      % true = fftshift (zero Doppler in the middle)
P.fft.chunk       = 64;         % FFT frames processed per block (memory control)

%% ------------------------------------------------------------------ Non-coherent integration
P.nci.nFrames = 4;              % n : buffer length (consecutive FFT outputs per bin)
P.nci.law     = 'square';       % detector law: 'square' |x|^2 | 'linear' |x| | 'log' dB
P.nci.average = true;           % true = mean over the buffer, false = sum
P.nci.output  = 'valid';        % 'valid' : only full buffers
                                % 'same'  : keep all frames (partial buffers at start)

%% ------------------------------------------------------------------ CFAR
P.cfar.type          = 'SO';    % 'SO' (smallest-of) | 'CA' | 'GO'
P.cfar.nRef          = 16;      % reference cells on EACH side
P.cfar.nGuard        = 3;       % guard cells on EACH side
P.cfar.thresholdMode = 'pfa';   % 'pfa'    : multiplier computed from P.cfar.pfa
                                % 'factor' : multiplier given by P.cfar.factorDb
P.cfar.pfa           = 1e-6;    % design probability of false alarm
P.cfar.factorDb      = 13;      % threshold factor [dB] when thresholdMode = 'factor'
P.cfar.nIntEffective = [];      % effective looks per integrated cell (Pfa design)
                                % [] = computed from the canceler, FFT window,
                                %      hop and nci.nFrames (see rsp_cfar_looks)
P.cfar.nRefEffective = [];      % effective independent cells per reference window
                                % [] = computed from the matched filter, whose
                                %      output noise is correlated in range
P.cfar.edge          = 'oneSided'; % range edges where one window is incomplete:
                                % 'oneSided' : use the complete side only
                                % 'partial'  : use the available cells
                                % 'none'     : no detection there
P.cfar.bins          = [];      % Doppler bins to test ([] = all)
P.cfar.rangeCells    = [];      % [first last] range cells to test ([] = all)
end
