%MAIN_IF Improvement factor of the MTI, measured on recorded clutter.
%
%   A few hundred consecutive pulses of raw video are enough. You do not need
%   to know where the clutter is: rsp_improvement_factor finds the clutter
%   cells (a stationary echo well above the noise) and the noise cells (no
%   pulse-to-pulse correlation) by itself.
%
%   Steps
%     1. input   : workspace (main_mcps.m), a .mat file, or simulated clutter
%     2. params  : radar_params.m (the two pulses, PRF, rpm, canceler, FFT)
%     3. measure : rsp_improvement_factor -> printed summary and figure
%
%   I = CNR before the canceler / CNR after it, the IEEE improvement factor
%   (S/C gain averaged over all target velocities). It is also given per range
%   cell, per pulse zone (short pulse near the radar, long pulse beyond), for
%   2/3/4-pulse cancelers and for every Doppler bin of the FFT that follows.
clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

%% 1. Input ---------------------------------------------------------------------------
% 'workspace' : "video" in the workspace, pulses x range, complex
%               (e.g. run main_mcps.m first)
% 'mat'       : a .mat file with one complex pulses x range matrix
%               (separate I and Q matrices: video = complex(I, Q))
% 'simulate'  : 500 simulated pulses, clutter up to 12 km and a target in it
source   = 'workspace';
matFile  = '';           % 'mat': file ('' = pick with a dialog)
matVar   = '';           % 'mat': variable ('' = first complex matrix in the file)
rowsAreRange = false;    % true if the rows of your matrix are range cells
dataType = 'video';      % 'video'   : raw video, decoded here with both pulses
                         % 'decoded' : already the decoder output (pulse compressed)
                         % 'samples' : raw samples, not decoded (codes not known)

P = radar_params();

switch source
    case 'workspace'
        if ~exist('video', 'var')
            error('main_if:input', 'No "video" in the workspace: run main_mcps.m first.');
        end
    case 'mat'
        if isempty(matFile)
            [f, d] = uigetfile('*.mat', 'Raw video (pulses x range, complex)');
            if isequal(f, 0)
                return
            end
            matFile = fullfile(d, f);
        end
        s = load(matFile);
        if isempty(matVar)
            names = fieldnames(s);
            for i = 1:numel(names)
                x = s.(names{i});
                if isnumeric(x) && ~isreal(x) && ismatrix(x) && min(size(x)) > 1
                    matVar = names{i};
                    break
                end
            end
            if isempty(matVar)
                error('main_if:input', 'No complex matrix in %s.', matFile);
            end
        end
        video = s.(matVar);
        clear s x
        fprintf('%s: "%s", %d x %d\n', matFile, matVar, size(video));
    case 'simulate'
        S = struct('nRange', 2500, 'noiseDb', 0, 'blankTx', true, 'seed', 1);
        S.targets = struct('rangeM', 8000, 'azDeg', 4, 'velocityMps', 20, 'powerDb', 10);
        S.clutter = struct('powerDb', 60, 'maxRangeM', 12000, 'sigmaVMps', 2, 'textureDb', 0);
        video = rsp_simulate(P, S, (1:500)');
        % expected: the simulated clutter has a Gaussian spectrum of std sigmaVMps
        g   = rsp_geometry(P, S.nRange);
        sf  = 2 * S.clutter.sigmaVMps / g.lambdaM / P.radar.prfHz;    % cycles/pulse
        fprintf('Simulated clutter, spread %.1f m/s. Expected I:', S.clutter.sigmaVMps);
        for order = [2 3 4]
            C   = P.canceler;
            C.order = order;
            w   = rsp_canceler_taps(C);
            rho = exp(-2*pi^2 * sf^2 * (0:numel(w)-1).^2);
            fprintf('  %d-pulse %.1f dB', order, 10*log10(sum(w.^2) / (w * toeplitz(rho) * w')));
        end
        fprintf('\n');
    otherwise
        error('main_if:input', 'Unknown source "%s".', source);
end
if rowsAreRange
    video = video.';
end
if size(video, 1) > size(video, 2)
    warning('main_if:layout', ['%d rows x %d columns: the rows should be pulses. ' ...
            'Set rowsAreRange = true if they are range cells.'], size(video));
end

%% 2.-3. Improvement factor -------------------------------------------------------------
IF = rsp_improvement_factor(video, P, 'input', dataType);

% Options, e.g.:
%   'cells', [70 800]          only these range cells (e.g. the clutter area)
%   'noiseCells', [3000 5000]  a range known to hold noise only
%   'noise', 10^(40/10)        the noise power per cell, if you know it
%   'minCnrDb', 20             only strong clutter
%   'maxVelocityMps', 2        only clutter close to zero Doppler
%   'pulses', 1:300            part of the record
% IF = rsp_improvement_factor(video, P, 'input', dataType, 'cells', [70 800]);
