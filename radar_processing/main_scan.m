%MAIN_SCAN Simulated rotating radar: full scan, live PPI, plots against truth.
%
%   The scenario is given in physical units: target range [m], azimuth [deg],
%   radial velocity [m/s] (positive = approaching) and SNR. The video is
%   generated pulse by pulse as the antenna turns (PRF, rpm, antenna
%   pattern), processed by the same chain as the real data (radar_params.m),
%   and shown on a PPI with a rotating sweep.
clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

P = radar_params();

%% Scenario ------------------------------------------------------------------------------
S.nRange     = 5469;            % range cells per pulse (fs = 6 MHz -> 25 m, 137 km)
S.noisePower = 1;
S.blankTx    = true;            % receiver off while transmitting
S.seed       = 1;

% SNR = per sample at the beam peak, before compression. The decoders add
% 10*log10(code length * samples per chip) (Barker 13: 19 dB, MLS 63: 26 dB)
% and the Doppler FFT about 10 dB.
S.targets = struct( ...
    'rangeM',      { 6000,  18000, 40000,  65000,  90000, 110000, 125000}, ...
    'azDeg',       {   30,     75,   120,    200,    250,    310,    340}, ...
    'velocityMps', {   40,    -20,   -75,    150,    200,     12,    -95}, ...
    'snrDb',       {  -10,    -12,   -15,    -16,    -18,    -15,    -20});

% Weak ground clutter, all azimuths, up to 25 km
S.clutter = struct('cnrDb', 20, 'maxRangeM', 25000, 'sigmaVMps', 0.5, 'textureDb', 3);

%% Run ------------------------------------------------------------------------------------
% P.output.bins = 2:16;          % e.g. leave out zero Doppler from the output
% P.ppi.source  = 'decoder';     % raw view instead of the output after CFAR
scan = rsp_scan(P, S, 'degrees', 360, 'display', true);

%% Replay ---------------------------------------------------------------------------------
% Real antenna speed (6 rpm = 10 s per turn); 'bins' shows chosen bins only:
% rsp_ppi_replay(scan, 'speed', 1);
% rsp_ppi_replay(scan, 'speed', 4, 'bins', [6 7 8]);
