clear;
clc;

file = '';
lanes = {};
scan = 1;
offset = 0;
count = 10000;

% read_mcps(file, lanes, options...)
%   file        ''                       pick with a dialog; any part (_p2, _p3...) finds the others
%   lanes       {}                       all of {'video','decoder','canceler','integral','cfar'}
%   'offset'    0                        PRFs to skip in the first loaded lane
%   'count'     Inf                      PRFs to read after the offset (Inf = everything)
%   'align'     true                     cut the other lanes to the same time span as the first lane
%   'crc'       true                     check every record's CRC-32 (false = faster)
%   'progress'  true                     show the progress bar
%   'quiet'     false                    true = do not print the summary table
%   'maxAzm'    8192                     azimuth counts per turn, used for scan numbering without marks
log = read_mcps(file, lanes, 'offset', offset, 'count', count);

% log.summary()                          print the per-lane table again
% log.details(lanes, options...)         gaps, out-of-order, lost records, lanes missing a pulse, field mismatches
%   'show'      20                       rows printed per table
%   'quiet'     false                    true = return the tables without printing
details = log.details();

% log.plot(scan, lanes, options...)      flat view: range cell x azimuth
% log.polar(scan, lanes, options...)     polar view, north up
%   'doppler'       'max'                integral / cfar: 'max' over bins, or a bin number (1..16)
%   'db'            true                 false = linear magnitude
%   'clim'          []                   fixed colour range, e.g. [40 110]
%   'dynamicRange'  []                   show the strongest value and this many dB below it, e.g. 50
%   'low'           50                   automatic contrast: lower percentile (noise floor)
%   'high'          99.95                automatic contrast: upper percentile
%   'range'         []                   range cells to show, e.g. [70 2500]
%   'colormap'      'jet'                any MATLAB colormap: 'parula', 'gray', 'hot', 'turbo'...
%   'maxPulses'     2000                 thin the rows to at most this many pulses
%   'maxCells'      1500                 thin the columns, keeping each block's maximum
log.polar(scan);
log.plot(scan);
log.plot(scan, {'decoder', 'canceler'});
log.plot(scan, {'decoder'}, 'doppler', 'max');
% log.plot(scan, {'decoder'}, 'dynamicRange', 50, 'range', [70 2500]);
% log.polar(scan, {'cfar'}, 'doppler', 5, 'colormap', 'hot');
% log.plot(scan, {'canceler'}, 'clim', [40 110], 'colormap', 'gray');

% log.matrix(lane, scans, channel)       scans [] = all loaded; channel for raw video (default 1)
%   video / decoder / canceler           pulses x range, complex (i + jq)
%   integral / cfar                      pulses x range x doppler, single
%   info: azimuth, azimuthDeg, timestampUs, pulseSeq, prfNumber, mode, scan, rangeCell (+ dopplerBin)
[video, videoInfo] = log.matrix('video', scan);
[decoder, decoderInfo] = log.matrix('decoder', scan);
[canceler, cancelerInfo] = log.matrix('canceler', scan);
[integral, integralInfo] = log.matrix('integral', scan);
[cfar, cfarInfo] = log.matrix('cfar', scan);

% log.export(file, lanes, scans)         saves each lane's matrix and info to a -v7.3 .mat file
[folder, name] = fileparts(log.file);
log.export(fullfile(folder, sprintf('%s_scan%d.mat', name, scan)), {}, scan);
