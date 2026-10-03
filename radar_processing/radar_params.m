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
P.radar.rpm        = 6;             % antenna speed [rev/min]
P.radar.fcMHz      = 1300;          % carrier [MHz] (only for m/s <-> Doppler)
P.radar.azStartDeg = 0;             % azimuth of the first pulse [deg]

P.antenna.beamwidthDeg = 4;         % 3 dB azimuth beamwidth [deg]
P.antenna.sidelobeDb   = -25;       % sidelobes, one way [dB]
P.antenna.backlobeDb   = -40;       % back lobe, one way [dB]

%% Sampling -----------------------------------------------------------------------------
P.fs = 6;                           % [MHz]

%% Pulse 1: short pulse ---------------------------------------------------------------------
% Phase code: put your chip values in .code and your decoder taps in .decoder
% (leave .decoder = [] for the matched filter). One chip = 1/bwMHz us.
P.pulse(1).type        = 'code';                     % 'code' | 'lfm'
P.pulse(1).bwMHz       = 1;                          % 1 MHz -> 1 us chips
P.pulse(1).code        = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];   % Barker 13
P.pulse(1).decoder     = [];                         % e.g. your mismatched filter
P.pulse(1).decoderForm = 'fir';                      % 'fir' (convolution taps) | 'reference'
P.pulse(1).decoderRate = 'chip';                     % 'chip' | 'sample'
P.pulse(1).fcMHz       = -1.5;                       % 0 if both pulses share one frequency
P.pulse(1).delayUs     = 63;                         % starts when pulse 2 ends
P.pulse(1).widthUs     = 13;                         % used by 'lfm' only

%% Pulse 2: long pulse ----------------------------------------------------------------------
P.pulse(2).type        = 'code';
P.pulse(2).bwMHz       = 1;
P.pulse(2).code        = rsp_code('mls', 63);        % 63-chip m-sequence
P.pulse(2).decoder     = [];
P.pulse(2).decoderForm = 'fir';
P.pulse(2).decoderRate = 'chip';
P.pulse(2).fcMHz       = 1.5;
P.pulse(2).delayUs     = 0;                          % sent first
P.pulse(2).widthUs     = 63;                         % used by 'lfm' only
P.pulse(2).window      = 'none';                     % e.g. 'taylor' for an LFM

% LFM instead of a code (1 MHz, Taylor weighted):
%   P.pulse(2).type = 'lfm';  P.pulse(2).widthUs = 63;
%   P.pulse(2).window = 'taylor';  P.pulse(2).windowParam = [4 -35];

%% Processing -----------------------------------------------------------------------------
P.canceler.order = 3;               % 3-pulse canceler
P.fft.nPulses    = 16;
P.fft.nfft       = 16;
P.fft.hop        = 1;
P.nci.nFrames    = 4;               % n
P.cfar.type      = 'SO';
P.cfar.nRef      = 48;
P.cfar.nGuard    = 8;
P.cfar.pfa       = 1e-6;

%% Output and display -------------------------------------------------------------------------
P.output.bins    = [];              % [] = max over all bins, e.g. 5 or [2:16]
P.ppi.source     = 'max';           % 'max' | 'video' | 'decoder' | 'canceler'
P.ppi.ringKm     = 20;
end
