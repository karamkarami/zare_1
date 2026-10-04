%MAIN_CHAIN Processing chain of one MCPS scan, compared with the log lanes.
%
%   video -> decoder -> 3-pulse canceler -> Doppler FFT -> NCI buffer -> SO-CFAR -> max
%
%   1. data       run main_mcps.m first (video, videoInfo, decoder, canceler,
%                 integral, cfar ... in the workspace), or load the .mat file of
%                 log.export, or make a synthetic test signal
%   2. parameters mc_params.m (put the decoder coefficients there)
%   3. chain      every block, one after the other
%   4. compare    a) my chain from the video against every log lane
%                 b) every block fed with the previous LOG lane (stage-wise)
%   5. figures
clc;
addpath(fileparts(mfilename('fullpath')));

%% 1. Data --------------------------------------------------------------------------------
source  = 'workspace';      % 'workspace' : run main_mcps.m first
                            % 'mat'       : the file written by log.export (matFile)
                            % 'synthetic' : test signal, no log needed
matFile = '';               % '' = choose in a dialog

%% 2. Parameters --------------------------------------------------------------------------
P = mc_params();

ref = struct('decoder', [], 'decoderInfo', [], 'canceler', [], 'cancelerInfo', [], ...
             'integral', [], 'integralInfo', [], 'cfar', [], 'cfarInfo', []);
truth = [];
switch source
    case 'workspace'
        if ~exist('video', 'var')
            error('main_chain:data', 'No "video" in the workspace: run main_mcps.m first.');
        end
        if ~exist('videoInfo', 'var'),    videoInfo = [];                end
        if exist('decoder', 'var'),       ref.decoder      = decoder;      end
        if exist('decoderInfo', 'var'),   ref.decoderInfo  = decoderInfo;  end
        if exist('canceler', 'var'),      ref.canceler     = canceler;     end
        if exist('cancelerInfo', 'var'),  ref.cancelerInfo = cancelerInfo; end
        if exist('integral', 'var'),      ref.integral     = integral;     end
        if exist('integralInfo', 'var'),  ref.integralInfo = integralInfo; end
        if exist('cfar', 'var'),          ref.cfar         = cfar;         end
        if exist('cfarInfo', 'var'),      ref.cfarInfo     = cfarInfo;     end
    case 'mat'
        [video, videoInfo, ref] = mc_load_mat(matFile);
    case 'synthetic'
        [video, videoInfo, P, truth] = mc_synthetic(P);
    otherwise
        error('main_chain:data', 'Unknown source "%s".', source);
end

fprintf('Video: %d pulses x %d range cells, %s\n', size(video, 1), size(video, 2), class(video));
az = mc_info(videoInfo, 'azimuthDeg');
if ~isempty(az)
    fprintf('  azimuth %.2f .. %.2f deg\n', min(az), max(az));
end
modes = unique(mc_info(videoInfo, 'mode'));
if ~isempty(modes)
    fprintf('  modes: %s\n', mat2str(modes(:).'));
end
lanes = {'decoder', 'canceler', 'integral', 'cfar'};
for i = 1:numel(lanes)
    if ~isempty(ref.(lanes{i}))
        fprintf('  log %-8s : %s\n', lanes{i}, mat2str(size(ref.(lanes{i}))));
    end
end

%% 3. Chain -------------------------------------------------------------------------------
rows0 = mc_rows(videoInfo, size(video, 1));          % row bookkeeping (pulseSeq, valid)

% 3.1 decoder: each pulse with its mismatched filter, short + long stitched at switchCell
tic;
[dec, decEach, filt] = mc_decoder(video, P);
elapsed.decoder = toc;
fprintf('\nDecoder (switch at cell %d):\n', P.decoder.switchCell);
for k = 1:numel(filt)
    fprintf('  pulse %d %-5s %4d samples, %4d taps, peak lag %4d (%s), shift %4d', k, ...
            P.pulse(k).name, filt(k).nSamples, filt(k).nCoef, filt(k).lag, ...
            filt(k).lagFrom, filt(k).shift);
    if ~isnan(filt(k).pslDb)
        fprintf(', PSL %.1f dB, loss %.2f dB', filt(k).pslDb, -filt(k).lossDb);
    end
    if filt(k).placeholder
        fprintf('   << PLACEHOLDER: put the coefficients in mc_params.m');
    end
    fprintf('\n');
end

% 3.2 3-pulse canceler
tic;
[can, rowsC] = mc_canceler(dec, P, rows0);
elapsed.canceler = toc;

% 3.3 Doppler FFT (on the same pulses as the log integral rows, when they are loaded)
tic;
endPulses = mc_fft_ends(rowsC, ref.integralInfo, size(ref.integral, 1), P);
[dop, rowsD] = mc_doppler(can, P, rowsC, endPulses);
elapsed.doppler = toc;
if isempty(endPulses)
    fprintf('Doppler FFT: one FFT every %d pulses -> %d rows\n', P.doppler.hop, size(dop, 1));
else
    fprintf('Doppler FFT: on the pulses of the %d log integral rows -> %d rows\n', ...
            size(ref.integral, 1), size(dop, 1));
end

% 3.4 non-coherent integration (buffer of P.nci.length FFT outputs on every bin)
tic;
[nci, rowsN] = mc_nci(dop, P, rowsD);
clear dop
elapsed.nci = toc;

% 3.5 SO-CFAR along range on every Doppler bin
tic;
cfarOut = mc_cfar(nci, P);
elapsed.cfar = toc;

% 3.6 max over the Doppler bins
[maxOut, maxBin] = mc_max(cfarOut, P);
fprintf('Time [s]: decoder %.1f, canceler %.1f, FFT %.1f, NCI %.1f, CFAR %.1f\n', ...
        elapsed.decoder, elapsed.canceler, elapsed.doppler, elapsed.nci, elapsed.cfar);
fprintf('Detections after the max: %d cells in %d rows\n', nnz(maxOut), size(maxOut, 1));

if ~isempty(truth)                                   % synthetic: were the targets found?
    for i = 1:numel(truth)
        c   = truth(i).cell + (-1:1);
        hit = maxOut(:, c);
        [~, j] = max(hit(:));
        b   = maxBin(:, c);
        fprintf('  target cell %4d bin %2d: detected in %3d of %d rows, bin of the max %d\n', ...
                truth(i).cell, truth(i).bin, nnz(any(hit > 0, 2)), size(maxOut, 1), b(j));
    end
end

%% 4. Comparison with the log --------------------------------------------------------------
Pq = P;
Pq.compare.plot = false;                             % numbers only
if ~isempty(ref.decoder) || ~isempty(ref.canceler) || ~isempty(ref.integral) || ~isempty(ref.cfar)
    fprintf('\na) My chain from the video against the log\n');
    mc_compare();
    mc_compare(dec, rows0, ref.decoder, ref.decoderInfo, 'decoder', P);
    mc_compare(can, rowsC, ref.canceler, ref.cancelerInfo, 'canceler', P);
    Rint = mc_compare(nci, rowsN, ref.integral, ref.integralInfo, 'integral', P);
    mc_compare(cfarOut, rowsN, ref.cfar, ref.cfarInfo, 'cfar', Pq, Rint.shift);
    if ~isempty(ref.cfar)
        mc_compare(maxOut, rowsN, max(ref.cfar, [], 3), ref.cfarInfo, 'max (final output)', ...
                   P, Rint.shift);
    end
    if ~isempty(ref.decoder)
        mc_align_check(decEach, rows0, ref.decoder, ref.decoderInfo, P, filt);
    end

    if P.compare.stageWise
        fprintf('\nb) Each block fed with the previous log lane\n');
        mc_compare();
        if ~isempty(ref.decoder) && ~isempty(ref.canceler)
            r = mc_rows(ref.decoderInfo, size(ref.decoder, 1));
            [x, r] = mc_canceler(ref.decoder, P, r);
            mc_compare(x, r, ref.canceler, ref.cancelerInfo, 'canceler <- log decoder', Pq);
        end
        if ~isempty(ref.canceler) && ~isempty(ref.integral)
            r = mc_rows(ref.cancelerInfo, size(ref.canceler, 1));
            e = mc_fft_ends(r, ref.integralInfo, size(ref.integral, 1), P);
            [x, r] = mc_doppler(ref.canceler, P, r, e);
            [x, r] = mc_nci(x, P, r);
            mc_compare(x, r, ref.integral, ref.integralInfo, 'integral <- log canceler', Pq);
        end
        if ~isempty(ref.integral) && ~isempty(ref.cfar)
            r = mc_rows(ref.integralInfo, size(ref.integral, 1));
            x = mc_cfar(ref.integral, P);
            mc_compare(x, r, ref.cfar, ref.cfarInfo, 'cfar <- log integral', Pq, 0);
        end
        clear x r e
    end
end
clear decEach

%% 5. Figures --------------------------------------------------------------------------------
figure('Name', 'My chain', 'NumberTitle', 'off');
subplot(2, 2, 1); mc_show(dec, 'decoder');
subplot(2, 2, 2); mc_show(can, 'canceler');
subplot(2, 2, 3); mc_show(nci, 'integral (max over bins)');
subplot(2, 2, 4); mc_show(maxOut, 'after SO-CFAR, max over bins');
