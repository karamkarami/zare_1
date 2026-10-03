%MAIN_RSP Run the processing chain on MCPS video and check every block.
%
%   video -> matched filter -> 3-pulse canceler -> Doppler FFT
%         -> non-coherent integration -> SO-CFAR
%
%   Steps
%     1. input   : the MCPS log (read_mcps), the workspace of main_mcps.m,
%                  or synthetic data
%     2. params  : every value of the chain, edit them here
%     3. chain   : full chain from the video
%     4. compare : a) full chain against every log lane
%                  b) every block fed with the previous LOG lane, so a
%                     mismatch points at one block only
%     5. plots
clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

%% 1. Input ---------------------------------------------------------------------------
% 'workspace' : use video / decoder / canceler / integral / cfar already in
%               the workspace (run main_mcps.m first)
% 'log'       : read the log here (needs read_mcps on the path)
% 'simulate'  : synthetic targets + clutter + noise (no reference lanes)
source = 'workspace';
scan   = 1;

switch source
    case 'workspace'
        if ~exist('video', 'var')
            error('main_rsp:input', 'No "video" in the workspace: run main_mcps.m first.');
        end
    case 'log'
        mcps = read_mcps('', {}, 'offset', 0, 'count', 10000);
        [video, videoInfo] = mcps.matrix('video', scan);
        decoder  = mcps.matrix('decoder', scan);
        canceler = mcps.matrix('canceler', scan);
        integral = mcps.matrix('integral', scan);
        cfar     = mcps.matrix('cfar', scan);
    case 'simulate'
        decoder = []; canceler = []; integral = []; cfar = [];
    otherwise
        error('main_rsp:input', 'Unknown source "%s".', source);
end
% lanes missing from the workspace are skipped in the comparison
if ~exist('decoder', 'var'),  decoder  = []; end
if ~exist('canceler', 'var'), canceler = []; end
if ~exist('integral', 'var'), integral = []; end
if ~exist('cfar', 'var'),     cfar     = []; end

%% 2. Parameters ----------------------------------------------------------------------
P = rsp_default_params();

% Sampling and the two pulses (p1, p2 in us, fs in MHz)
P.fs = 4;

P.pulse(1).widthUs     = 2;          % p1 : short pulse
P.pulse(1).bwMHz       = 0;          % 0 = unmodulated
P.pulse(1).fcMHz       = -1.2;
P.pulse(1).delayUs     = 20;         % starts after the long pulse
P.pulse(1).window      = 'none';

P.pulse(2).widthUs     = 20;         % p2 : long LFM pulse
P.pulse(2).bwMHz       = 2;
P.pulse(2).fcMHz       = 0.8;
P.pulse(2).slope       = +1;
P.pulse(2).delayUs     = 0;
P.pulse(2).window      = 'taylor';
P.pulse(2).windowParam = [4 -35];

% Matched filter
P.mf.combine    = 'stitch';          % 'stitch' | 'pulse1' | 'pulse2'
P.mf.switchCell = [];                % [] = blind zone of the long pulse
P.mf.norm       = 'noise';

% Canceler
P.canceler.order  = 3;
P.canceler.output = 'same';

% Doppler FFT
P.fft.nPulses = 16;
P.fft.nfft    = 16;
P.fft.hop     = 1;
P.fft.window  = 'hamming';

% Non-coherent integration
P.nci.nFrames = 4;                   % n consecutive outputs per bin
P.nci.law     = 'square';
P.nci.average = true;

% CFAR
P.cfar.type          = 'SO';
P.cfar.nRef          = 16;
P.cfar.nGuard        = 3;
P.cfar.thresholdMode = 'pfa';        % or 'factor' with P.cfar.factorDb
P.cfar.pfa           = 1e-6;
P.cfar.factorDb      = 13;

if strcmp(source, 'simulate')
    S.nPulses    = 350;
    S.nRange     = 5469;
    S.noisePower = 1;
    S.blankTx    = true;
    S.seed       = 1;
    S.targets    = struct('rangeCell', {50, 1000, 2500, 4000}, ...
                          'fdNorm',    {0.25, 3/16, -0.31, 0.02}, ...
                          'snrDb',     {0, -10, -15, -10});
    S.clutter    = struct('cnrDb', 40, 'cells', [1 1500], 'spreadNorm', 0.01);
    video = rsp_simulate(P, S);
end

%% 3. Full chain ----------------------------------------------------------------------
out = rsp_chain(video, P);
fprintf('Chain: %d x %d video, switch cell %d, %d detections\n', size(video), ...
        out.mf.switchCell, numel(out.cfar.list.range));
fprintf('Time [s]: decoder %.2f  canceler %.2f  doppler %.2f  integral %.2f  cfar %.2f\n', ...
        out.time.decoder, out.time.canceler, out.time.doppler, out.time.integral, out.time.cfar);
fprintf('CFAR factor %s dB, effective looks %s, effective ref cells %s\n\n', ...
        mat2str(10*log10(out.cfar.factor), 3), mat2str(out.nInt, 3), mat2str(out.nRef, 3));

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
        blkCfar = rsp_cfar(integral, P.cfar, P.nci.law, out.nInt, out.nRef);
        rsp_compare(cfar, blkCfar.map, 'cfar      <- log integral');
    end

    % Visual check of one lane, e.g.:
    % rsp_compare(decoder, out.decoder, 'decoder', 'plot', true);
end

%% 5. Plots ----------------------------------------------------------------------------
rsp_plot(out, video);
% rsp_plot(out, video, 'range', [1 1500], 'bin', 5);
