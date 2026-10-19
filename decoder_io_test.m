clc
clear
close all

%% ============================ PARAMETERS ============================
fs = 2e6;                
c  = 3e8;             
noise_on = true;
rng(1);

% ----------------------------- fixed point -----------------------------
% decoder input  : int16 (complex), one ADC scale for the whole PRI
% decoder taps   : int16, final FIR taps (window / flip / conj already applied), scaled to full range
% decoder output : int32 = (full-precision accumulator) >> out_shift, out_shift from worst-case growth
adc_peak = [];            % input amplitude mapped to int16 full scale ([] = peak of rx in the PRI)

% ----------------------------- short pulse -----------------------------
sp.type = 'phase';        % 'phase' or 'lfm'
sp.chip_us  = 1;          % chip width [us] -> round(chip_us*1e-6*fs) samples per chip
sp.code     = [1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 -1 1 -1 -1 1 1 1 1 1];
mmf = [-0.170847936550266;
0.006480974129943;
0.021239501875478;
0.021312912832743;
0.116597689016265;
-0.115023766372816;
-0.012992001717144;
-0.032781669281148;
0.073105097047959;
0.261451725938703;
-0.311516909948656;
0.193502924275457;
-0.429405800686943;
0.355393283430697;
-0.071144083306238;
0.222472375040466;
-0.186380471992882;
-0.223384692725169;
-0.197564011383335
0.412598793095109;
-0.213467723468599;
1.134344563012496;
-1.632186792867563;
0.811369838456468;
-1.264640688107512;
1.385644074052266;
0.671010106609803;
-0.765960061276992;
-0.407353133419193;
-0.887457583593748;
0.755037617760453;
1.569991081687012;
-0.755037617760452;
-0.887457583593748;
0.407353133419193;
-0.765960061276992;
-0.671010106609803;
1.385644074052267;
1.264640688107512;
0.811369838456468;
1.632186792867563;
1.134344563012496;
0.213467723468599;
0.412598793095109;
0.197564011383336;
-0.223384692725169;
0.186380471992882;
0.222472375040466;
0.071144083306238;
0.355393283430697;
0.429405800686943;
0.193502924275457;
0.311516909948656;
0.261451725938703;
-0.073105097047959;
-0.032781669281148;
0.012992001717144;
-0.115023766372816;
-0.116597689016265;
0.021312912832743;
-0.021239501875479;
0.006480974129943;
0.170847936550266];
sp.mismatch = mmf;        % one tap per chip (3 x code length); repeated to samples per chip internally
sp.mismatch_form = 'correlator';   % 'correlator': taps in code order (y = xcorr(x, taps))
                                   % 'impulse'   : taps already are the FIR impulse response (y = conv(x, taps))
sp.lfm_T  = 100e-6;        % pulse width [s]
sp.lfm_B  = 1e6;          % sweep bandwidth [Hz]
sp.window = 'none';     % 'taylor', 'hamming' or 'none'
sp.taylor_nbar = 4;       % Taylor: number of nearly constant sidelobes
sp.taylor_sll  = -35;     % Taylor: sidelobe level [dB]
sp.gap_us    = 10;         % gap1: end of short TX -> short listening starts [us]
sp.listen_us = 60;       % short listening window [us]

sp.target_range_m    = 5000;    
sp.target_snr_db     = 30;     
sp.target_doppler_hz = 0;     

gap_mid_us = 30;           % gap2: end of short listening -> long TX starts [us]

lp.type = 'lfm';        % 'phase' or 'lfm'
lp.chip_us  = 1;
lp.code     = kron([1 1 1 -1 1], [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]);
lp.mismatch = [];
lp.mismatch_form = 'correlator';
lp.lfm_T  = 100e-6;
lp.lfm_B  = 1e6;
lp.window = 'taylor';
lp.taylor_nbar = 4;
lp.taylor_sll  = -35;
lp.gap_us    = 10;         % gap3: end of long TX -> long listening starts [us]
lp.listen_us = 250;       % long listening window [us]
lp.target_range_m    = 30000;  
lp.target_snr_db     = 30;
lp.target_doppler_hz = 0;

sp.name = 'Short';
lp.name = 'Long';
[sp.tx, sp.h, sp.filter_delay] = make_pulse_and_filter(sp, fs);
[lp.tx, lp.h, lp.filter_delay] = make_pulse_and_filter(lp, fs);

sp.tx_start  = 0;
sp.win_start = sp.tx_start + numel(sp.tx) + round(sp.gap_us*1e-6*fs);
sp.win_len   = round(sp.listen_us*1e-6*fs);
lp.tx_start  = sp.win_start + sp.win_len + round(gap_mid_us*1e-6*fs);
lp.win_start = lp.tx_start + numel(lp.tx) + round(lp.gap_us*1e-6*fs);
lp.win_len   = round(lp.listen_us*1e-6*fs);
N_pri        = lp.win_start + lp.win_len;

fprintf('PRI = %d samples (%.1f us): short TX %d | gap %d | short listen %d | gap %d | long TX %d | gap %d | long listen %d\n\n', ...
        N_pri, N_pri/fs*1e6, numel(sp.tx), sp.win_start - numel(sp.tx), sp.win_len, ...
        lp.tx_start - sp.win_start - sp.win_len, numel(lp.tx), lp.win_start - lp.tx_start - numel(lp.tx), lp.win_len);

sp.echo_start = sp.tx_start + round(2*sp.target_range_m/c*fs);
lp.echo_start = lp.tx_start + round(2*lp.target_range_m/c*fs);

rx = zeros(1, N_pri);
if noise_on
    rx = (randn(1, N_pri) + 1j*randn(1, N_pri)) / sqrt(2);
end
rx = add_echo(rx, sp, fs, c);
rx = add_echo(rx, lp, fs, c);

listening = false(1, N_pri);
listening(sp.win_start + (1:sp.win_len)) = true;
listening(lp.win_start + (1:lp.win_len)) = true;
rx(~listening) = 0;                       % receiver is off outside the listening windows

% ---------------- fixed point: int16 input, int16 taps, int32 output ----------------
if isempty(adc_peak), adc_peak = max(abs([real(rx), imag(rx)])); end
in_scale = double(intmax('int16')) / adc_peak;           % rx units -> int16 counts
[sp.hq, sp.h_scale, sp.out_shift] = quantize_filter(sp.h);
[lp.hq, lp.h_scale, lp.out_shift] = quantize_filter(lp.h);

dec_in_short  = int16(rx(sp.win_start + (1:sp.win_len)) * in_scale);
dec_in_long   = int16(rx(lp.win_start + (1:lp.win_len)) * in_scale);
dec_out_short = decoder(dec_in_short, sp.hq, sp.out_shift);
dec_out_long  = decoder(dec_in_long,  lp.hq, lp.out_shift);

h_short = sp.hq;  h_scale_short = sp.h_scale;  out_shift_short = sp.out_shift;  delay_short = sp.filter_delay;
h_long  = lp.hq;  h_scale_long  = lp.h_scale;  out_shift_long  = lp.out_shift;  delay_long  = lp.filter_delay;
tx_short = sp.tx;  tx_long = lp.tx;
save('decoder_io.mat', 'dec_in_short', 'dec_out_short', 'dec_in_long', 'dec_out_long', ...
     'h_short', 'h_long', 'h_scale_short', 'h_scale_long', 'out_shift_short', 'out_shift_long', ...
     'delay_short', 'delay_long', 'in_scale', 'tx_short', 'tx_long', 'fs');

fprintf('Fixed point: input int16 (scale %.2f counts/unit) | taps int16 | output int32\n', in_scale);
fprintf('   short: %d taps, tap scale %.1f, output shift %d\n', numel(sp.hq), sp.h_scale, sp.out_shift);
fprintf('   long : %d taps, tap scale %.1f, output shift %d\n\n', numel(lp.hq), lp.h_scale, lp.out_shift);

% back to rx units (noise power 1) for checks and plots
x_short = double(dec_in_short) / in_scale;
x_long  = double(dec_in_long)  / in_scale;
y_short = double(dec_out_short) * 2^sp.out_shift / (in_scale * sp.h_scale);
y_long  = double(dec_out_long)  * 2^lp.out_shift / (in_scale * lp.h_scale);

%% =================== 5) Check peaks and plot ===================
check_peak(y_short, sp, fs, c, noise_on);
check_peak(y_long,  lp, fs, c, noise_on);
plot_results(rx, sp, lp, x_short, y_short, x_long, y_long, fs, c);


%% ============================ FUNCTIONS ============================
function y = decoder(x, h, out_shift)
% DECODER BLOCK under test: FIR filter y(n) = sum_k h(k) x(n-k+1) (full convolution).
% x: int16 input, h: int16 final taps, y: int32 = accumulator >> out_shift (rounded, saturated).
% length(y) = length(x) + length(h) - 1. The peak of an echo that starts at
% input sample k appears at output sample k + filter_delay.
    acc = conv(double(x), double(h));                  % exact: all values are integers < 2^53
    y   = int32(acc / 2^out_shift);
end

function [hq, scale, out_shift] = quantize_filter(h)
% Final taps -> int16 (largest real/imag part = 32767), and the output shift that
% keeps the worst-case accumulator (|input| = 32768 on every tap) inside int32.
    scale = double(intmax('int16')) / max(abs([real(h), imag(h)]));
    hq    = int16(h * scale);
    hd    = double(hq);
    acc_max   = 32768 * sum(abs(real(hd)) + abs(imag(hd)));
    out_shift = max(0, ceil(log2(acc_max / double(intmax('int32')))));
end

function [tx, h, delay] = make_pulse_and_filter(p, fs)
% Builds the transmitted baseband pulse and the final decoder taps h
% (window, flip and conj already applied: the decoder only does conv(x, h)).
    switch p.type
        case 'phase'
            spc  = round(p.chip_us*1e-6*fs);                    % samples per chip
            code = p.code(:).';
            tx   = repelem(code, spc);
            if isempty(p.mismatch)
                h_chip = ls_mismatch_filter(code, 3*numel(code));   % placeholder until you insert yours
            else
                h_chip = p.mismatch(:).';
                if strcmp(p.mismatch_form, 'correlator')
                    h_chip = conj(fliplr(h_chip));              % correlator taps -> impulse response
                end
            end
            if numel(h_chip) ~= 3*numel(code)
                warning('%s: mismatch filter has %d taps, expected 3*%d = %d.', ...
                        p.name, numel(h_chip), numel(code), 3*numel(code));
            end
            if max(abs(conv(code, conj(fliplr(h_chip))))) > max(abs(conv(code, h_chip)))
                warning('%s: mismatch taps look time-reversed. Check mismatch_form.', p.name);
            end
            h = repelem(h_chip, spc);                           % spc taps per chip (6x code length at 2 samples/chip)
        case 'lfm'
            N  = round(p.lfm_T*fs);
            t  = ((0:N-1) - (N-1)/2) / fs;                       % time centred on the pulse
            tx = exp(1j*pi*(p.lfm_B/p.lfm_T)*t.^2);              % up-chirp, -B/2 .. +B/2
            h  = conj(fliplr(tx)) .* make_window(p, N);          % matched filter x window
        otherwise
            error('%s: type must be ''phase'' or ''lfm''.', p.name);
    end
    tx = tx / sqrt(mean(abs(tx).^2));                            % unit average power (scale only)

    % Decoder delay: where the peak appears for an echo that starts at sample 1
    [~, ipk] = max(abs(conv(tx, h)));
    delay = ipk - 1;
end

function w = make_window(p, N)
% Amplitude window for the LFM matched filter.
    switch p.window
        case 'taylor',  w = taylorwin(N, p.taylor_nbar, p.taylor_sll).';
        case 'hamming', w = hamming(N).';
        otherwise,      w = ones(1, N);
    end
end

function h = ls_mismatch_filter(code, M)
% Least-squares mismatch filter of length M (only a placeholder design).
    N = numel(code);
    A = convmtx(code(:), M);              % conv(code, h) = A*h
    d = zeros(N+M-1, 1);
    d(ceil((N+M-1)/2)) = 1;               % desired response: one spike in the middle
    h = (A \ d).';
end

function rx = add_echo(rx, p, fs, c)
% Adds the target echo (delayed, Doppler-shifted copy of the pulse) to rx.
    n    = 0:numel(p.tx)-1;
    amp  = 10^(p.target_snr_db/20);                              % noise power is 1
    echo = amp * p.tx .* exp(1j*2*pi*p.target_doppler_hz*n/fs);
    if p.echo_start < p.win_start || p.echo_start + numel(p.tx) > p.win_start + p.win_len
        r_min = (p.win_start - p.tx_start) / fs * c/2;
        r_max = (p.win_start + p.win_len - numel(p.tx) - p.tx_start) / fs * c/2;
        warning('%s target should be between %.0f m and %.0f m (fully inside its window).', ...
                p.name, r_min, r_max);
    end
    idx = p.echo_start + n + 1;                                  % 1-based indices in rx
    ok  = idx <= numel(rx);
    rx(idx(ok)) = rx(idx(ok)) + echo(ok);
end

function check_peak(y, p, fs, c, noise_on)
% Prints where the decoder output peak is and how good the decoder filter is.
    [pk, k]  = max(abs(y));
    k_expect = (p.echo_start - p.win_start + 1) + p.filter_delay;
    range_m  = (p.win_start + k - 1 - p.filter_delay - p.tx_start) / fs * c/2;
    [pslr_db, loss_db] = filter_quality(p.tx, p.h);

    fprintf('%s pulse (%s): pulse %d samples, decoder %d taps, PSLR %.1f dB, mismatch loss %.2f dB\n', ...
            p.name, p.type, numel(p.tx), numel(p.h), pslr_db, loss_db);
    fprintf('   peak at output sample %d (expected %d) -> range %.0f m (target %.0f m)', ...
            k, k_expect, range_m, p.target_range_m);
    if abs(k - k_expect) <= 1, fprintf('   OK\n'); else, fprintf('   WRONG\n'); end
    if noise_on
        snr_out    = 20*log10(pk / norm(p.h));                   % output noise power = sum|h|^2
        snr_theory = p.target_snr_db + 10*log10(numel(p.tx)) + loss_db;
        fprintf('   output peak SNR %.1f dB (theory %.1f dB)\n', snr_out, snr_theory);
    end
    fprintf('\n');
end

function [pslr_db, loss_db] = filter_quality(tx, h)
% Noise-free decoder response to the clean pulse.
% PSLR = highest sidelobe / peak. Loss = SNR loss versus the ideal matched filter.
    r = abs(conv(tx, h));
    [pk, ipk] = max(r);
    a = ipk;  while a > 1        && r(a-1) < r(a), a = a - 1; end   % main lobe left edge
    b = ipk;  while b < numel(r) && r(b+1) < r(b), b = b + 1; end   % main lobe right edge
    sidelobes = r([1:a-1, b+1:end]);
    pslr_db = 20*log10(max([sidelobes, eps]) / pk);
    loss_db = 20*log10(pk / (norm(tx)*norm(h)));
end

function plot_results(rx, sp, lp, in_s, out_s, in_l, out_l, fs, c)
    figure('Name', 'Decoder test', 'Position', [80 60 1250 850]);

    % --- PRI timeline ---
    subplot(3, 2, [1 2]); hold on; grid on;
    t_us = (0:numel(rx)-1) / fs * 1e6;
    tx_line = zeros(size(rx));
    tx_line(sp.tx_start + (1:numel(sp.tx))) = abs(sp.tx);
    tx_line(lp.tx_start + (1:numel(lp.tx))) = abs(lp.tx);
    win_a = [sp.win_start, lp.win_start] / fs * 1e6;
    win_b = [sp.win_start + sp.win_len, lp.win_start + lp.win_len] / fs * 1e6;
    xregion(win_a, win_b, 'FaceColor', [0.3 0.6 1]);
    h1 = plot(t_us, abs(rx), 'Color', [0.45 0.45 0.45]);
    h2 = plot(t_us, tx_line, 'r', 'LineWidth', 1.5);
    top = 1.1 * max([abs(rx), 1]);
    text((win_a + win_b)/2, [top top], {'short listen', 'long listen'}, 'HorizontalAlignment', 'center');
    ylim([0 1.2*top]);
    legend([h1 h2], {'RX (shaded = listening)', 'TX pulses'}, 'Location', 'northeast');
    xlabel('time in PRI [\mus]'); ylabel('|signal|');
    title('PRI: short TX | gap1 | short listen | gap2 | long TX | gap3 | long listen');

    % --- decoder input / output ---
    plot_io(3, in_s, out_s, sp, fs, c);
    plot_io(5, in_l, out_l, lp, fs, c);
end

function plot_io(pos, x, y, p, fs, c)
% Decoder input (left) and output (right) versus range.
    to_km = @(n) (n - p.tx_start) / fs * c/2 / 1e3;   % sample -> range of an echo starting there
    r_in  = to_km(p.win_start + (0:numel(x)-1));
    r_out = to_km(p.win_start + (0:numel(y)-1) - p.filter_delay);

    subplot(3, 2, pos);
    plot(r_in, abs(x)); grid on;
    xline(p.target_range_m/1e3, 'r--');
    title(sprintf('%s window: decoder INPUT (%s)', p.name, p.type));
    xlabel('range [km]'); ylabel('|x|');

    subplot(3, 2, pos + 1);
    plot(r_out, 20*log10(abs(y)/max(abs(y)) + eps)); grid on;
    xline(p.target_range_m/1e3, 'r--');
    ylim([-70 3]);
    title(sprintf('%s window: decoder OUTPUT', p.name));
    xlabel('range [km]'); ylabel('dB (re peak)');
end
