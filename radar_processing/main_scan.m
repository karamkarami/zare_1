%MAIN_SCAN Simulated rotating radar: live PPI, plots against truth, real-time replay.
%
%   The scenario is given in physical units: target range [m], azimuth [deg],
%   radial velocity [m/s] (positive = approaching) and echo power [dB]. The
%   video is generated pulse by pulse as the antenna turns (PRF, rpm,
%   antenna pattern), processed by the same chain as the real data
%   (radar_params.m), and shown on a PPI whose sweep turns with the antenna.
%
%   Steps
%     1. link budget : expected SNR of each target through every block and
%                      its margin over the CFAR threshold (rsp_budget)
%     2. live scan   : simulate + process; the PPI sweep paints as it goes
%     3. report      : every target echo against its plot (errors, margin)
%     4. replay      : the processed revolutions at the real antenna speed
clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

P = radar_params();

%% Scenario ------------------------------------------------------------------------------
% All levels are powers per sample at the receiver, in dB on one scale.
S.nRange  = 5469;               % range cells per pulse (fs = 6 MHz -> 25 m, 137 km)
S.noiseDb = 0;                  % receiver noise
S.blankTx = true;               % receiver off while transmitting
S.seed    = 1;

%                 range [m]  az [deg]  velocity [m/s]  power [dB]
targets = [         10000       40          20            10
                    50000       20          60             5
                    80000       80         100             0 ];
S.targets = struct('rangeM',      num2cell(targets(:, 1))', ...
                   'azDeg',       num2cell(targets(:, 2))', ...
                   'velocityMps', num2cell(targets(:, 3))', ...
                   'powerDb',     num2cell(targets(:, 4))');

% Ground clutter, all azimuths, up to 25 km
S.clutter = struct('powerDb', -50, 'maxRangeM', 25000, 'sigmaVMps', 0.5, 'textureDb', 3);

nTurns = 2;                     % antenna revolutions (10 s each at 6 rpm)

%% Run ------------------------------------------------------------------------------------
% P.output.bins = 2:16;          % e.g. leave zero Doppler out of the output
% P.ppi.source  = 'decoder';     % raw view instead of the output after CFAR
% P.ppi.speed   = Inf;           % do not wait for the real antenna speed
scan = rsp_scan(P, S, 'degrees', 360 * nTurns, 'display', true);

%% Replay at the real antenna speed (6 rpm = 10 s per turn) ---------------------------------
rsp_ppi_replay(scan, 'speed', 1, 'loops', 1);
% rsp_ppi_replay(scan, 'speed', 4, 'bins', [4 5 6]);    % faster, chosen bins only
