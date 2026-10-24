%MAIN_PPI 2-D rotating search radar (PPI) in terrain clutter, with a scan schedule.
%
%   A rotating surveillance radar is simulated pulse by pulse (no tracking):
%   noise, ground / mountain / forest clutter, terrain shadow and moving
%   targets. Every scan of the schedule uses its own carrier frequency, PRF,
%   2- or 3-pulse canceler, Doppler FFT length (pulses integrated coherently)
%   or output bins, and they follow each other on one live PPI.
%
%   Steps
%     1. scenario  : clutter map, targets and antenna pattern (rsp_plot_scenario)
%     2. schedule  : the scans, live PPI (rsp_schedule -> rsp_scan)
%     3. numbers   : per scan the link budget, every target against its plot
%                    (errors in azimuth / range / velocity), false alarms per
%                    clutter region against the design Pfa; the timeline
%     4. figures   : summary of every scan (rsp_plot_summary)
%     5. save      : results/run_<date>/  .mat + CSV + text + PNG (rsp_save)
%   Run rsp_selftest to verify every block (validation tests).
clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

P = radar_params();             % carrier, PRF, rpm, waveform, canceler, FFT, CFAR ...

%% Antenna pattern --------------------------------------------------------------------------
% Built-in patterns: P.antenna.type = 'model' (default) | 'sinc' | 'gaussian'.
% Your own pattern, two columns angle [deg], one-way gain [dB] (.csv / .txt / .mat):
% P.antenna = rsp_antenna_load(P.antenna, fullfile(here, 'patterns', 'example_pattern.csv'));

%% Scenario -----------------------------------------------------------------------------------
% All levels are powers per sample at the receiver, in dB on one scale (noise = 0 dB).
S.nRange  = 5469;               % range cells (25 m each at fs = 6 MHz -> 137 km)
S.noiseDb = 0;                  % receiver noise
S.blankTx = true;               % receiver off while transmitting
S.seed    = 1;

%                 range [m]  az [deg]  velocity [m/s]  power [dB]
targets = [         10000       40          20            10      % 1 short-pulse zone, slow
                    50000       20          60             5      % 2
                    80000       80         100             0      % 3
                    35000      120          40             0      % 4 in front of the mountain
                    75000      125          80             5      % 5 behind the mountain (shadow)
                    20000      225         -50             0      % 6 over the forest, receding
                   105000      300         150            10 ];   % 7 fast: folds in velocity
S.targets = struct('rangeM',      num2cell(targets(:, 1))', ...
                   'azDeg',       num2cell(targets(:, 2))', ...
                   'velocityMps', num2cell(targets(:, 3))', ...
                   'powerDb',     num2cell(targets(:, 4))');

% Clutter layers (rsp_clutter_types documents every field and the defaults of each type)
S.clutter = struct('type', {'ground', 'mountain', 'forest'});
S.clutter(1).rangeM     = [0 20000];       % ground around the radar, all azimuths
S.clutter(1).powerDb    = 15;
S.clutter(2).azDeg      = [100 140];       % mountain: strong, almost static, hides what is behind
S.clutter(2).rangeM     = [40000 55000];
S.clutter(2).powerDb    = 45;
S.clutter(2).shadowDb   = -40;
S.clutter(3).azDeg      = [200 250];       % forest: weaker, moves in the wind
S.clutter(3).rangeM     = [12000 30000];
S.clutter(3).powerDb    = 25;
S.clutter(3).sigmaVMps  = 0.6;

%% Scan schedule ------------------------------------------------------------------------------
% Every entry starts from P and changes only what it lists; time and azimuth
% run on from one entry to the next (targets keep moving).
turns = 1;                      % revolutions per entry (10 s each at 6 rpm)
sched = struct( ...
  'name',  {'Search 600 MHz PRF 1000', '2-pulse canceler', 'FFT 32 pulses', ...
            'PRF 1500 Hz, 1300 MHz', 'Low Doppler bins off'}, ...
  'turns', {turns, turns, turns, turns, turns}, ...
  'set',   {{}, ...
            {'canceler.order', 2}, ...
            {'fft.nPulses', 32, 'fft.nfft', 32, 'nci.nFrames', 16}, ...
            {'radar.prfHz', 1500, 'radar.fcMHz', 1300}, ...
            {'output.bins', 3:15}});

% P.ppi.speed = Inf;            % do not wait for the real antenna speed
% P.ppi.source = 'decoder';     % raw video on the PPI instead of the CFAR output

%% Run ----------------------------------------------------------------------------------------
rsp_plot_scenario(P, S);
run = rsp_schedule(P, S, sched, 'display', true);

%% Figures and saved results -----------------------------------------------------------------
for k = 1:numel(run.scans)
    rsp_plot_summary(run.scans(k));
end
files = rsp_save(run, '');      % results/run_<date>/: results.mat, *.csv, summary.txt, *.png

% Replay one scan of the schedule at 4x the antenna speed:
%   rsp_ppi_replay(run.scans(1), 'speed', 4);
% Load the saved results again:
%   s = load(fullfile(fileparts(files{1}), 'results.mat'));  s.result.timeline
