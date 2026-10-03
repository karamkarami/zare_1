function P = radar_params()
%RADAR_PARAMS Values of YOUR radar. Edit this file; main_rsp and main_scan use it.
%
%   P = radar_params()
%
%   Starts from rsp_default_params (which documents every field) and sets
%   the values that describe the radar. Anything not set here keeps its
%   default.

P = rsp_default_params();

%% Radar --------------------------------------------------------------------------------
P.radar.prfHz      = 1000;          % PRF [Hz]
P.radar.rpm        = 6;             % antenna speed [rev/min] -> 10 s per turn
P.radar.fcMHz      = 600;           % carrier [MHz]; sets m/s <-> Doppler:
                                    % unambiguous velocity = +-lambda*PRF/4
                                    % (600 MHz: +-125 m/s; 1300 MHz: +-58 m/s)
P.radar.azStartDeg = 0;             % azimuth of the first pulse [deg]

P.antenna.beamwidthDeg = 4;         % 3 dB azimuth beamwidth [deg]
P.antenna.sidelobeDb   = -25;       % sidelobes, one way [dB]
P.antenna.backlobeDb   = -40;       % back lobe, one way [dB]

%% Sampling -----------------------------------------------------------------------------
P.fs = 6;                           % [MHz] -> 25 m range cells

%% Pulse 1: short pulse, Barker 13 ------------------------------------------------------
% Phase code: chip values in .code, your decoder taps in .decoder
% ([] = matched filter). One chip = 1/bwMHz us.
P.pulse(1).type        = 'code';                     % 'code' | 'lfm'
P.pulse(1).bwMHz       = 1;                          % 1 MHz -> 1 us chips
P.pulse(1).code        = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
P.pulse(1).decoder     = rsp_code_mmf(P.pulse(1).code, 52);  % mismatched, PSL -48 dB
P.pulse(1).decoderForm = 'fir';                      % 'fir' (convolution taps) | 'reference'
P.pulse(1).decoderRate = 'chip';                     % 'chip' | 'sample'
P.pulse(1).fcMHz       = -1.5;                       % 0 if both pulses share one frequency
P.pulse(1).delayUs     = 100;                        % starts when pulse 2 ends

%% Pulse 2: long pulse, LFM 100 us -------------------------------------------------------
P.pulse(2).type        = 'lfm';
P.pulse(2).widthUs     = 100;                        % [us]
P.pulse(2).bwMHz       = 1;                          % swept bandwidth [MHz]
P.pulse(2).slope       = 1;                          % +1 up, -1 down
P.pulse(2).fcMHz       = 1.5;
P.pulse(2).delayUs     = 0;                          % sent first
P.pulse(2).window      = 'taylor';                   % range sidelobes (PSL -41 dB)
P.pulse(2).windowParam = [6 -45];

% Phase code instead of the LFM:
%   P.pulse(2).type = 'code';  P.pulse(2).code = <your chips>;  P.pulse(2).window = 'none';
%   P.pulse(2).decoder = <your taps>;  P.pulse(1).delayUs = numel(P.pulse(2).code);

%% Processing -----------------------------------------------------------------------------
P.canceler.order = 3;               % 3-pulse canceler
P.fft.nPulses    = 16;
P.fft.nfft       = 16;
P.fft.hop        = 1;
P.nci.nFrames    = 32;              % n: 16 + 32 - 1 = 47 pulses (1.7 deg) < beam dwell
P.cfar.type      = 'SO';
P.cfar.nRef      = [];              % [] = 16 resolution cells (96 cells at 6 MHz)
P.cfar.nGuard    = [];              % [] = widest main lobe + 2 (14 cells here)
P.cfar.pfa       = 1e-6;

%% Output and display -------------------------------------------------------------------------
P.output.bins    = [];              % [] = max over all bins, e.g. 5 or [2:16]
P.ppi.source     = 'max';           % 'max' | 'video' | 'decoder' | 'canceler'
P.ppi.ringKm     = 20;
P.ppi.stepDeg    = 0.5;             % sweep step per screen update [deg]
P.ppi.speed      = 1;               % 1 = real antenna speed (6 rpm), Inf = fastest
end
