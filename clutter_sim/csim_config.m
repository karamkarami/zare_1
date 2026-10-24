function S = csim_config()
%CSIM_CONFIG Settings of the test-data simulator.
%
%   S = csim_config()
%
%   The simulator is independent of the processor: it writes complex baseband
%   after pulse compression (pulses x range), with the same block fields as
%   online data (x, prfHz, azDeg, rangeM, velEN), so the processor cannot
%   tell the two apart.
%
%   Echo of a point at range R(t): a * g(az) * exp(-j 4 pi R(t) / lambda),
%   so the Doppler of the clutter and of the targets comes from the geometry
%   (platform and target motion), g is the two-way Gaussian beam.

%% Radar
S.fcHz         = 600e6;
S.beamwidthDeg = 4;
S.rpm          = 6;
S.azStartDeg   = 0;
S.rangeCellM   = 150;          % range cell r is at r * rangeCellM
S.nRange       = 300;          % 45 km
S.noisePower   = 1;            % per sample
S.blockPulses  = 256;          % pulses per block (the last block of a scan stops at north)
S.seed         = 1;

%% Ground clutter: one scatterer per range cell every azStepDeg, fixed on the ground
S.clutter.maxRangeM  = 20e3;
S.clutter.minRangeM  = 600;
S.clutter.azStepDeg  = 0.5;
S.clutter.cnrDb      = 45;     % mean CNR per sample at refRangeM
S.clutter.refRangeM  = 5e3;
S.clutter.rangeExp   = 2;      % mean power ~ (R / refRangeM)^-rangeExp
S.clutter.textureDb  = 5;      % log-normal spread of the scatterer powers (std) [dB]
S.clutter.sigmaVMps  = 0.3;    % intrinsic spectrum spread (sum of components with Gaussian Doppler)
S.clutter.components = 4;
% patches of stronger clutter: [azimuth width(deg) range(m) depth(m) extra(dB)]
S.clutter.patches = [120 25 12e3 4e3 15;
                     250 30  6e3 3e3 10];

%% Platform (ground vehicle): velocity [east north] in m/s, can be changed between blocks
S.platform.velEN = [0 0];

%% Targets: start position, radial / tangential speed over ground, SNR per sample on boresight
%  vRadMps > 0 = approaching.
S.swerling = true;             % true: Swerling 1, amplitude drawn again every revolution
                               % false: constant amplitude (to compare with the predicted gains)
S.targets = struct( ...
    'rangeM',  {15e3,  8e3,   40e3,  10e3,  18e3,  4e3}, ...
    'azDeg',   {30,    120,   200,   300,   250,   80}, ...
    'vRadMps', {100,  -220,   300,   40,   -25,   -330}, ...
    'vTanMps', {0,     0,     0,     0,     0,     0}, ...
    'snrDb',   {15,    15,    15,    15,    15,    15});
end
