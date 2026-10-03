%MAIN_RSP Process one block of MCPS video and check every block of the chain.
%
%   video -> decoder -> 3-pulse canceler -> Doppler FFT -> integration
%         -> SO-CFAR -> max over bins -> plots
%
%   Steps
%     1. input   : workspace of main_mcps.m, the MCPS log, or simulated data
%     2. params  : radar_params.m (your codes, decoders, fs, PRF ...)
%     3. chain   : full chain from the video
%     4. compare : a) full chain against every log lane
%                  b) every block fed with the previous LOG lane, so a
%                     mismatch points at one block only
%     5. output  : plots (range m, azimuth deg, velocity m/s), figures, PPI sector
clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

%% 1. Input ---------------------------------------------------------------------------
% 'workspace' : video / decoder / canceler / integral / cfar already in the
%               workspace (run main_mcps.m first)
% 'log'       : read the log here (needs read_mcps on the path)
% 'simulate'  : 350 simulated pulses (targets, weak clutter, antenna pattern)
source  = 'workspace';
scanNum = 1;

P = radar_params();

switch source
    case 'workspace'
        if ~exist('video', 'var')
            error('main_rsp:input', 'No "video" in the workspace: run main_mcps.m first.');
        end
    case 'log'
        mcps = read_mcps('', {}, 'offset', 0, 'count', 10000);
        [video, videoInfo] = mcps.matrix('video', scanNum);
        decoder  = mcps.matrix('decoder', scanNum);
        canceler = mcps.matrix('canceler', scanNum);
        integral = mcps.matrix('integral', scanNum);
        cfar     = mcps.matrix('cfar', scanNum);
    case 'simulate'
        % 350 pulses (12.6 deg) around the three targets, moved next to each
        % other in azimuth so that one block sees all of them
        S = struct('nRange', 5469, 'noiseDb', 0, 'blankTx', true, 'seed', 1);
        S.targets = struct('rangeM',      {10000, 50000, 80000}, ...
                           'azDeg',       {4,     6,     8}, ...
                           'velocityMps', {20,    60,    100}, ...
                           'powerDb',     {10,    5,     0});
        S.clutter = struct('powerDb', -50, 'maxRangeM', 25000, 'sigmaVMps', 0.5, 'textureDb', 3);
        video = rsp_simulate(P, S, (1:350)');
        rsp_budget(P, S);
        decoder = []; canceler = []; integral = []; cfar = [];
        tr = rsp_truth(P, S);
        fprintf('Simulated targets: az [deg]  range [m]  vel [m/s]  measured vel [m/s]  bin\n');
        for k = 1:numel(tr)
            fprintf('                  %8.2f %10.0f %10.1f %14.1f %8d\n', tr(k).azDeg, ...
                    tr(k).rangeM, tr(k).velocityMps, tr(k).foldedMps, tr(k).bin);
        end
    otherwise
        error('main_rsp:input', 'Unknown source "%s".', source);
end
% lanes missing from the workspace are skipped in the comparison
if ~exist('decoder', 'var'),  decoder  = []; end
if ~exist('canceler', 'var'), canceler = []; end
if ~exist('integral', 'var'), integral = []; end
if ~exist('cfar', 'var'),     cfar     = []; end

%% 3. Full chain ----------------------------------------------------------------------
out = rsp_chain(video, P);
fprintf('Chain: %d x %d video, switch cell %d, valid cells %d..%d\n', size(video), ...
        out.mf.switchCell, out.mf.validCells);
for k = 1:numel(P.pulse)
    d = out.mf.dec{k};
    fprintf('  pulse %d: %-4s %4d samples  PSL %6.1f dB  ISL %6.1f dB  loss %4.2f dB  lag %d\n', ...
            k, P.pulse(k).type, numel(d.tx), d.pslDb, d.islDb, abs(d.lossDb), d.lag);
    if d.lossDb < -3
        warning('main_rsp:decoder', ['Pulse %d decoder loses %.1f dB against a matched ' ...
                'filter: check decoderForm (''fir'' / ''reference'') and decoderRate.'], k, -d.lossDb);
    end
end
fprintf('Time [s]: decoder %.2f  canceler %.2f  doppler %.2f  integral %.2f  cfar %.2f  max %.2f\n', ...
        out.time.decoder, out.time.canceler, out.time.doppler, out.time.integral, ...
        out.time.cfar, out.time.max);
fprintf('CFAR %s: %d guard + %d reference cells, factor %s dB\n\n', out.P.cfar.type, ...
        out.P.cfar.nGuard, out.P.cfar.nRef, mat2str(10*log10(out.cfar.factor), 3));

%% 4. Comparison with the log lanes -----------------------------------------------------
if ~isempty(decoder) || ~isempty(canceler) || ~isempty(integral) || ~isempty(cfar)
    % Rows of the 3-D lanes: if the log keeps one row per pulse, take the
    % rows that match the pulse index of our frames.
    pick = @(lane, idx) lane(idx(idx <= size(lane, 1)), :, :);
    integralRef = integral;
    cfarRef     = cfar;
    if size(integral, 1) == size(video, 1)
        integralRef = pick(integral, out.idx.integral);
    end
    if size(cfar, 1) == size(video, 1)
        cfarRef = pick(cfar, out.idx.cfar);
    end

    fprintf('a) Full chain from video against the log\n');
    rsp_compare(decoder,     out.decoder,  'decoder');
    rsp_compare(canceler,    out.canceler, 'canceler');
    rsp_compare(integralRef, out.integral, 'integral');
    rsp_compare(cfarRef,     out.cfar.map, 'cfar');

    fprintf('\nb) Each block fed with the previous log lane\n');
    if ~isempty(decoder) && ~isempty(canceler)
        blkCanceler = rsp_canceler(decoder, P.canceler);
        rsp_compare(canceler, blkCanceler, 'canceler  <- log decoder');
    end
    if ~isempty(canceler) && ~isempty(integral)
        [blkDoppler, idxD]  = rsp_doppler_fft(canceler, P.fft);
        [blkIntegral, idxI] = rsp_nci(blkDoppler, P.nci, idxD);
        if size(integral, 1) == size(canceler, 1)
            integralRef = pick(integral, idxI);
        end
        rsp_compare(integralRef, blkIntegral, 'integral  <- log canceler');
        clear blkDoppler
    end
    if ~isempty(integral) && ~isempty(cfar)
        if P.cfar.validOnly
            blkCfar = rsp_cfar_cells(integral, out.P, out.mf.validCells, out.nInt, out.nRef);
        else
            blkCfar = rsp_cfar(integral, out.P.cfar, P.nci.law, out.nInt, out.nRef);
        end
        rsp_compare(cfar, blkCfar.map, 'cfar      <- log integral');
    end
    % Visual check of one lane, e.g.:
    % rsp_compare(decoder, out.decoder, 'decoder', 'plot', true);
end

%% 5. Output -----------------------------------------------------------------------------
% pulse numbers and azimuth of the rows (from the log when available)
g = rsp_geometry(P, size(video, 2), (1:size(video, 1))');
if exist('videoInfo', 'var') && isfield(videoInfo, 'azimuthDeg')
    g.azDeg = double(videoInfo.azimuthDeg(:));
end

% plots from the max output
[f, c] = find(out.max.value > 0);
lin    = sub2ind(size(out.max.value), f, c);
det    = struct('pulse', out.idx.max(f), 'cell', c, ...
                'bin', double(out.max.bin(lin)), 'value', double(out.max.value(lin)), ...
                'binFrac', double(out.max.binFrac(lin)), 'marginDb', double(out.max.marginDb(lin)), ...
                'azDeg', g.azDeg(out.idx.max(f)));
plots  = rsp_extract_plots(det, P, size(video, 2));
fprintf('\n%d plots\n   az [deg]   range [m]   vel [m/s]  bin  CFAR margin [dB]  hits\n', numel(plots.azDeg));
for i = 1:numel(plots.azDeg)
    fprintf('  %8.2f  %10.0f   %8.1f  %3d  %15.1f  %5d\n', plots.azDeg(i), plots.rangeM(i), ...
            plots.velocityMps(i), plots.bin(i), plots.marginDb(i), plots.hits(i));
end

rsp_plot(out, video);
% rsp_plot(out, video, 'bin', 7, 'bins', [2:16], 'range', [1 2000]);

% PPI of this block (one sector)
ppi  = rsp_ppi_init(P, g.rangeM, 'PPI - this block');
rows = double(out.max.marginDb);                        % dB over the CFAR threshold
rows(~isfinite(rows)) = -Inf;
ppi  = rsp_ppi_update(ppi, g.azDeg(out.idx.max), rows, plots, 'Output after CFAR + max');
