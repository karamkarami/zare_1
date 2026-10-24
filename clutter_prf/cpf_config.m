function cfg = cpf_config()
%CPF_CONFIG Default settings of the clutter-map MTI and the scan-to-scan PRF set.
%
%   cfg = cpf_config()   edit the returned struct (or this file) for your radar.
%
%   The processor works on complex baseband (I + jQ) after pulse compression,
%   one block of consecutive pulses at a time. A block may be any length, but
%   all its pulses must have the same PRF (the PRF changes only between scans):
%
%       blk.x       pulses x range, complex
%       blk.prfHz   PRF of these pulses [Hz]
%       blk.azDeg   antenna azimuth of every pulse [deg], 0 = north, clockwise
%       blk.rangeM  range of every range cell [m]
%       blk.velEN   platform velocity [vEast vNorth] [m/s], 1x2 or pulses x 2
%                   (optional, default [0 0] = stationary radar)
%       blk.posEN   platform position [east north] [m], 1x2 or pulses x 2
%                   (optional; integrated from velEN when missing)
%
%   Online data and the simulator (clutter_sim) fill the same struct.

%% Radar
cfg.radar.fcHz         = 600e6;     % carrier frequency [Hz]; only sets lambda
cfg.radar.beamwidthDeg = 4;         % azimuth 3 dB beamwidth [deg]
cfg.radar.rpm          = 6;         % antenna speed, used until it is measured from blk.azDeg

%% Targets
cfg.target.vMinMps = 20;            % radial speed band over ground, both signs:
cfg.target.vMaxMps = 350;           % +-[vMin vMax] (approaching and receding)

%% PRF set (one PRF per scan, used in turn)
cfg.prf.minHz         = 800;
cfg.prf.maxHz         = 1500;
cfg.prf.stepHz        = 5;          % candidate grid (PRF resolution of the radar)
cfg.prf.count         = 4;          % P
cfg.prf.dvMps         = 0.5;        % velocity grid of the criteria [m/s]
cfg.prf.blindMode     = 'clutter';  % blind = 'clutter': folded Doppler inside the clutter region
                                    %         (clutter edge + half a Doppler FFT bin)
                                    %         'notch': filter gain below blindDb
cfg.prf.fftLength     = 16;         % Doppler FFT after the filter (sets the clutter region)
cfg.prf.blindDb       = -10;        % 'notch': gain re noise gain [dB]
cfg.prf.maxBlind      = 1;          % at any velocity at most this many of the P PRFs may be blind
cfg.prf.caMinDb       = [];         % clutter attenuation needed; [] = design CNR (clutter down to the noise)
cfg.prf.wRipple       = 1;          % cost per dB of ripple of the mean response (mean over the P PRFs)
cfg.prf.wCa           = 1;          % cost per dB the worst PRF's clutter attenuation misses caMinDb
cfg.prf.wPd           = 0;          % cost per % of min cumulative Pd over P scans (0 = only reported)
cfg.prf.snrDb         = 15;         % target SNR at 0 dB filter gain, for the Pd figures
cfg.prf.pfa           = 1e-6;
cfg.prf.maxExhaustive = 5e4;        % try every subset up to this many, otherwise local search
cfg.prf.restarts      = 10;         % local search starts

%% MTI filter (replaces the fixed canceler; the Doppler FFT follows it)
cfg.mti.maxLength   = 8;            % maximum number of taps
cfg.mti.lengthMode  = 'max';        % 'max' | 'auto' = shortest length whose residue <= residueDb
cfg.mti.residueDb   = -3;           % 'auto': clutter left after the filter, re noise [dB]
cfg.mti.cnrMarginDb = 6;            % design for this much more clutter than the map holds
cfg.mti.widthMargin = 1.2;          % design for a spectrum this much wider than the map holds
cfg.mti.designIter  = 150;          % iterations of the passband shaping (mean dB gain)
cfg.mti.bypassCnrDb = 0;            % cells with less clutter pass unfiltered
cfg.mti.vmStepMps   = 0.25;         % quantisation of the clutter mean velocity
cfg.mti.motionComp  = true;         % remove the platform Doppler of the ground (TACCAR)
cfg.mti.output      = 'ground';     % 'ground': platform Doppler stays removed | 'radar'
cfg.mti.sGrid       = logspace(-3.3, log10(0.3), 33);   % filter bank: clutter sigma / PRF (19% steps)
cfg.mti.cnrGridDb   = 0:3:90;                            % filter bank: CNR [dB]

%% Clutter map (fixed to the ground, polar around the map origin)
cfg.map.azCellDeg     = 2;          % about half the beamwidth
cfg.map.rangeCellM    = 600;        % a few range resolution cells
cfg.map.maxRangeM     = 60e3;
cfg.map.alpha         = 0.3;        % weight of the newest scan (memory of about 1/alpha scans)
cfg.map.clipDb        = 10;         % largest change of a cell in one update after the first [dB]
cfg.map.minUpdates    = 2;          % updates before a cell gets its own filter
cfg.map.minCnrDbWidth = 15;         % spectrum width and mean are measured only above this CNR
cfg.map.maxLag        = 32;         % correlation lags used for the spectrum width
cfg.map.recenterM     = 3000;       % moving platform: re-centre the map after this distance

%% Clutter assumed where the map has no data yet (and for the first PRF set)
cfg.clutter.cnrDb        = 60;
cfg.clutter.sigmaVMps    = 0.5;     % intrinsic spectrum spread (wind-blown vegetation)
cfg.clutter.sigmaVMinMps = 0.02;

%% Noise
cfg.noise.power    = [];            % noise power per sample; [] = estimated from the data
cfg.noise.fraction = 0.2;           % estimate: mean power of the range cells below twice the
                                    % power of the 20% weakest one (needs 20% of cells without clutter)
cfg.noise.alpha    = 0.1;

%% Periodic re-selection of the PRF set (cpf_next_prf)
cfg.update.everyScans      = 10;
cfg.update.speedDeltaMps   = 3;     % also re-select when the platform speed changed this much
cfg.update.hysteresis      = 0.5;   % switch only if the cost drops by this much
cfg.update.cnrPercentile   = 95;    % design clutter taken from the map
cfg.update.sigmaPercentile = 90;
cfg.update.motionCase      = 'worst';   % platform broadening: 'worst' (broadside) | 'mean'
end
