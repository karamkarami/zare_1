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



%% decoder
% video is pulses x range: each ROW is one pulse (fast time = range cells), each COLUMN is one
% range cell over the pulses (slow time).
% The decoder (pulse compression) works on every pulse along RANGE, i.e. along dim 2.
% conv2(A,h) with a ROW vector h convolves every row of A along dim 2, which is what we need; a
% COLUMN vector would convolve along the pulses. So the taps are forced to be a row vector with
% (:).' whatever shape they were saved in (.' is the plain transpose, ' would also conjugate them).
% decoder_up is the reference the echo is correlated with. Correlation with c is convolution with
% the time-reversed CONJUGATE of c, conj(flip(c)). For real taps (binary code) conj does nothing,
% for complex taps (polyphase code) it is needed.
% (If decoder_up is already the FIR impulse response and not the reference, drop flip and conj.)
load Decoder_35_1_up
mismatch_35 = conj(flip(decoder_up(:).'));
load Decoder_240_1_up
mismatch_240 = conj(flip(decoder_up240(:).'));

[n_pulse,n_range] = size(video);
n_short = 556;                                          % range cells 1..n_short: short pulse, the rest: long pulse

figure,mesh(real(video)),title('video')
video_short = video(:,1:n_short);
video_long = video(:,n_short+1:end);
decoder_out_short = conv2(video_short,mismatch_35);     % pulses x (n_short + length(mismatch_35) - 1)
decoder_out_long = conv2(video_long,mismatch_240);      % pulses x (n_range - n_short + length(mismatch_240) - 1)

% conv2 gives the 'full' convolution, longer than the input. An echo that starts at range cell r
% gives its compressed peak at cell r + delay:
%   matched filter (length L = code length)                    delay = L - 1
%   mismatch filter of length M for an L-sample code           delay = (L + M)/2 - 1  (peak in the middle)
% Dropping the first 'delay' cells puts the peak back on the cell where the echo starts and keeps
% the input width, so the short and long parts can be joined again.
% The check at the end of this file measures the shift against the radar's decoder lane; if it
% is not 0, add it to delay_35 / delay_240.
delay_35 = length(mismatch_35) - 1;
delay_240 = length(mismatch_240) - 1;
decoder_out_short = decoder_out_short(:,delay_35 + (1:n_short));
decoder_out_long = decoder_out_long(:,delay_240 + (1:n_range - n_short));
decoder_out = [decoder_out_short,decoder_out_long];     % pulses x range, same size as video
figure,mesh(20*log10(abs(decoder_out) + eps)),title('decoder [dB]')

%% canceller
% The 3-pulse canceller (MTI) works on every range cell along the PULSES (slow time), i.e. along
% dim 1. A row vector [-1/4 1/2 -1/4] in conv2 filters along range (dim 2), which is wrong: the
% taps must be a COLUMN vector, coef_canceller(:).
% Row n of the 'full' output is -1/4*x(n) + 1/2*x(n-1) - 1/4*x(n-2), the canceller at pulse n.
% Only rows 1..n_pulse are kept (the last 2 rows would need pulses after the last one), so every
% row stays on its own pulse, the same row as in video and videoInfo.
% Rows 1 and 2 are the start-up transient (the 3-pulse buffer is not full yet and the clutter is
% not cancelled), so they are set to 0.
% The filter is the same for every range cell, so the short and long parts are filtered together.
% The canceller and the FFT need consecutive pulses with the same PRF; if the scan mixes modes or
% PRFs (videoInfo.mode / videoInfo.prfNumber), process each group of pulses on its own.
coef_canceller = [-1/4 1/2 -1/4];
canceller_out = conv2(decoder_out,coef_canceller(:));
canceller_out = canceller_out(1:n_pulse,:);
canceller_out(1:length(coef_canceller)-1,:) = 0;
figure,mesh(20*log10(abs(canceller_out) + eps)),title('canceller [dB]')

%% FFT
% Sliding 16-pulse buffer on every range cell. At every pulse the newest canceller output goes
% in, the oldest drops out, and a 16-point FFT is taken along the pulses (dim 1 of the buffer).
% Buffer row 1 = oldest pulse (n-15), row 16 = newest pulse (n):
%   fft_out(n,r,k) = sum_{m=0}^{15} w(m+1) * canceller_out(n-15+m,r) * exp(-1j*2*pi*(k-1)*m/16)
% Doppler bin k = 1 is zero Doppler, bin k is (k-1)/16*PRF and bins 10..16 are the negative
% frequencies (k-17)/16*PRF (no fftshift: the same order 1..16 as the radar's lanes).
% The buffer starts with zeros, so the first 15 outputs hold fewer than 16 pulses.
n_fft = 16;
doppler_window = ones(n_fft,1);                         % rectangular = plain FFT; e.g. hamming(n_fft) for lower Doppler sidelobes
fft_buffer = complex(zeros(n_fft,n_range,class(canceller_out)));
fft_out = complex(zeros(n_pulse,n_range,n_fft,class(canceller_out)));
for n = 1:n_pulse
    fft_buffer = [fft_buffer(2:end,:); canceller_out(n,:)];      % oldest pulse out, newest pulse in
    spectrum = fft(fft_buffer.*doppler_window,n_fft,1);         % n_fft x n_range, FFT along the pulses
    fft_out(n,:,:) = reshape(spectrum.',1,n_range,n_fft);       % .' = plain transpose (no conjugate)
end

pulse_show = round(n_pulse/2);                          % pulse shown in the range-Doppler figures
figure,mesh(20*log10(abs(squeeze(fft_out(pulse_show,:,:))) + eps)),title(sprintf('FFT, pulse %d [dB]',pulse_show))
xlabel('Doppler bin'),ylabel('range cell')

%% non-coherent integration
% Every FFT output cell (range cell, Doppler bin) has its own 5-deep buffer with the last 5 FFT
% outputs. The phase is dropped (square-law detector |X|^2) and the 5 powers are summed:
%   integral_out(n,r,k) = sum_{i=0}^{4} |fft_out(n-i,r,k)|^2
% The 5 FFTs slide by one pulse, so one output covers 16 + 5 - 1 = 20 pulses.
% (For a linear detector use abs(...) instead of abs(...).^2.)
n_integ = 5;
integ_buffer = zeros(n_integ,n_range,n_fft,class(canceller_out));
integral_out = zeros(n_pulse,n_range,n_fft,class(canceller_out));
for n = 1:n_pulse
    integ_buffer = cat(1,integ_buffer(2:end,:,:),abs(fft_out(n,:,:)).^2);   % oldest FFT out, newest FFT in
    integral_out(n,:,:) = sum(integ_buffer,1);
end

% first pulse whose buffers hold only valid data: canceller (3) + FFT (16) + integration (5)
first_valid = length(coef_canceller) + n_fft + n_integ - 2;

figure,mesh(10*log10(max(integral_out,[],3) + eps)),title('integral, max over Doppler bins [dB]')
xlabel('range cell'),ylabel('pulse')
figure,mesh(10*log10(squeeze(integral_out(pulse_show,:,:)) + eps)),title(sprintf('integral, pulse %d [dB]',pulse_show))
xlabel('Doppler bin'),ylabel('range cell')

%% CFAR OS
% Ordered-Statistic CFAR along RANGE, on every pulse and on every FFT (Doppler) bin on its own:
% each bin is a separate range profile and its reference cells come from the same bin only.
% For the cell under test (CUT) r, with n_guard guard cells and n_ref reference cells on each side:
%   reference cells: r-n_guard-n_ref ... r-n_guard-1  and  r+n_guard+1 ... r+n_guard+n_ref
%   1. sort the reference cells (smallest first)
%   2. throw away the n_discard largest ones
%   3. noise level of the window = the largest cell that is left
%   4. threshold = alpha * noise level; detection if integral_out(r) > threshold
% Throwing away the largest cells keeps up to n_discard cells of other targets (or a clutter edge)
% in the window from raising the threshold.
% Near the ends of a part, only the reference cells inside the part are used (fewer cells), and
% alpha is taken for that number of cells, so the Pfa stays the same.
% After their decoders the short and long parts have different noise levels, so each part gets
% its own CFAR and no window crosses the switch cell (cfar_part = {1:n_range} for one CFAR).
% Pulses before first_valid hold the start-up transient of the buffers and are not tested.
n_guard = 2;                                            % guard cells on each side (must cover the main lobe of the compressed pulse)
n_ref = 16;                                             % reference cells on each side
n_discard = 5;                                          % largest reference cells thrown away
pfa = 1e-6;                                             % probability of false alarm
cfar_part = {1:n_short,n_short+1:n_range};

% alpha: the factor that turns the noise level into the threshold for the wanted Pfa.
% Pfa itself is not the multiplier: for N reference cells and the k-th smallest as the noise level
% (k = N - n_discard), square-law detector and Gaussian noise (Rohling, 1983):
%   Pfa = prod_{i=0}^{k-1} (N - i) / (N - i + alpha)
% This Pfa falls when alpha grows, so alpha is found by bisection, for every N that can occur.
% (The 5 integrated FFTs share 15 of their 16 pulses, so they are close to one look; this alpha
% is close to right and on the safe side, the real Pfa is a bit lower.)
alpha = nan(1,2*n_ref);                                 % alpha(N) for N reference cells
for N = n_discard+1:2*n_ref
    k = N - n_discard;
    i = 0:k-1;
    alpha_low = 0;                                      % Pfa = 1 here (above the wanted Pfa)
    alpha_high = 1e6;                                   % Pfa far below the wanted Pfa here
    for iter = 1:100
        alpha_mid = (alpha_low + alpha_high)/2;
        pfa_mid = prod((N - i)./(N - i + alpha_mid));
        if pfa_mid > pfa
            alpha_low = alpha_mid;                      % too many false alarms: raise alpha
        else
            alpha_high = alpha_mid;
        end
    end
    alpha(N) = (alpha_low + alpha_high)/2;
end
fprintf('CFAR OS: %d guard + %d reference cells per side, %d largest thrown away, Pfa %.0e\n',n_guard,n_ref,n_discard,pfa);
fprintf('         alpha = %.2f dB with %d reference cells (k = %d), %.2f dB with %d cells at the ends (k = %d)\n', ...
        10*log10(alpha(2*n_ref)),2*n_ref,2*n_ref-n_discard,10*log10(alpha(n_ref)),n_ref,n_ref-n_discard);

cfar_det = false(n_pulse,n_range,n_fft);
cfar_threshold = nan(n_pulse,n_range,n_fft,class(integral_out));
pulses = first_valid:n_pulse;
for p = 1:numel(cfar_part)
    cells = cfar_part{p};
    for r = cells
        % reference cells of this CUT (guard cells and the CUT itself are left out), inside the part
        ref_cells = [r-n_guard-n_ref:r-n_guard-1, r+n_guard+1:r+n_guard+n_ref];
        ref_cells = ref_cells(ref_cells >= cells(1) & ref_cells <= cells(end));
        N = numel(ref_cells);
        if N <= n_discard
            continue;                                   % too few reference cells: not tested
        end
        ref = sort(integral_out(pulses,ref_cells,:),2); % pulses x N x doppler, every pulse and bin sorted on its own
        noise = ref(:,N-n_discard,:);                   % largest cell left after the n_discard largest are thrown away
        threshold = alpha(N)*noise;                     % pulses x 1 x doppler
        cfar_threshold(pulses,r,:) = threshold;
        cfar_det(pulses,r,:) = integral_out(pulses,r,:) > threshold;
    end
end
cfar_out = integral_out.*cfar_det;                      % integrated value on detections, 0 elsewhere
fprintf('CFAR OS: %d detections (pulse x range x Doppler bin)\n',nnz(cfar_det));

[p_det,r_det] = find(any(cfar_det,3));
figure,plot(r_det,p_det,'.'),title('CFAR detections, any Doppler bin')
xlabel('range cell'),ylabel('pulse'),axis([1 n_range 1 n_pulse])

% CUT and threshold along range at the strongest detection (or the middle pulse, bin 1)
if any(cfar_det(:))
    [~,i_max] = max(cfar_out(:));
    [pulse_cfar,~,bin_cfar] = ind2sub(size(cfar_out),i_max);
else
    pulse_cfar = max(pulse_show,first_valid);
    bin_cfar = 1;
end
figure,plot(10*log10(integral_out(pulse_cfar,:,bin_cfar) + eps)),hold on
plot(10*log10(cfar_threshold(pulse_cfar,:,bin_cfar)),'r')
plot(find(cfar_det(pulse_cfar,:,bin_cfar)),10*log10(integral_out(pulse_cfar,cfar_det(pulse_cfar,:,bin_cfar),bin_cfar)),'ko')
legend('integral','threshold','detection'),xlabel('range cell'),ylabel('dB')
title(sprintf('CFAR OS, pulse %d, Doppler bin %d',pulse_cfar,bin_cfar))

%% compare with the radar's lanes
% The log holds the radar's own output of every block for the same scan. check_lane prints, for
% each block, the range shift that lines my output up with the radar's lane and how well they
% match there: the correlation of the log magnitudes (1 = same shape), so the radar's
% fixed-point scale and a |.| or |.|^2 detector do not change it.
% A positive shift means my output is that many cells late: add it to delay_35 / delay_240.
% A low correlation at the best shift points at the block itself (filter taps, direction ...).
if ~isempty(decoder)
    check_lane(decoder_out(:,1:n_short),decoder(:,1:min(n_short,end)),'decoder short');
    check_lane(decoder_out(:,n_short+1:end),decoder(:,n_short+1:end),'decoder long');
end
if ~isempty(canceler)
    check_lane(canceller_out,canceler,'canceller');
end
if ~isempty(integral)
    check_lane(max(integral_out,[],3),max(integral,[],3),'integral');   % max over bins: Doppler bin order does not matter
end
if ~isempty(cfar)
    n_p = min(n_pulse,size(cfar,1));
    n_r = min(n_range,size(cfar,2));
    radar_det = any(cfar(1:n_p,1:n_r,:) ~= 0,3);
    my_det = any(cfar_det(1:n_p,1:n_r,:),3);
    fprintf('%-14s radar %d, mine %d, both %d detected cells (pulse x range, any Doppler bin)\n','cfar',nnz(radar_det),nnz(my_det),nnz(radar_det & my_det));
end



%% local functions
function check_lane(mine,radar,name)
% Finds the range shift that lines mine up with the radar's lane and prints how well they match.
% mine, radar: pulses x range (complex or real). Only the common pulses / range cells are used,
% about 256 pulses spread over the scan to keep it fast.
max_shift = 64;                                         % range shifts tried: -max_shift ... +max_shift
n_p = min(size(mine,1),size(radar,1));
n_r = min(size(mine,2),size(radar,2));
rows = unique(round(linspace(1,n_p,min(n_p,256))));
a = log10(double(abs(mine(rows,1:n_r))));               % log magnitude: free of scale and of |.| / |.|^2
b = log10(double(abs(radar(rows,1:n_r))));
shift = -max_shift:max_shift;
rho = nan(size(shift));
for i = 1:numel(shift)
    r = max(1,1-shift(i)):min(n_r,n_r-shift(i));        % mine(:,r+shift) against radar(:,r)
    x = a(:,r+shift(i));
    y = b(:,r);
    ok = isfinite(x) & isfinite(y);                     % cells that are 0 (log = -Inf) are left out
    if nnz(ok) > 100
        c = corrcoef(x(ok),y(ok));
        rho(i) = c(1,2);
    end
end
[rho_best,i_best] = max(rho);
fprintf('%-14s best range shift %+3d cells: correlation %.3f (at shift 0: %.3f)\n',name,shift(i_best),rho_best,rho(shift == 0));
end
