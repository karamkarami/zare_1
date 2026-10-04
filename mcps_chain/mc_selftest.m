function ok = mc_selftest()
%MC_SELFTEST Check every block against a direct (slow) computation.
%
%   ok = mc_selftest()      prints PASS / FAIL for each check

nPass = 0;
nFail = 0;
    function check(name, cond)
        if cond
            nPass = nPass + 1;
            fprintf('  PASS  %s\n', name);
        else
            nFail = nFail + 1;
            fprintf('  FAIL  %s\n', name);
        end
    end

rng(7);
P0 = mc_params();
quiet = warning('off', 'mc_filter:placeholder');

% ---- decoder: peaks of point targets land on their range cells --------------------------------
S = struct('noiseDb', -300, 'clutterDb', -300, 'nPulses', 40, 'nCells', 2000, ...
           'targets', struct('cell', {300, 1200}, 'bin', {1, 1}, 'powerDb', {0, 0}));
[v, info, P] = mc_synthetic(P0, S);
[y, each] = mc_decoder(v, P);
a = abs(y(20, :));
[~, c1] = max(a(1:559));
[~, c2] = max(a(560:end));
check('decoder: short-pulse target at its cell', c1 == 300);
check('decoder: long-pulse target at its cell', c2 + 559 == 1200);
[~, c3] = max(abs(each{1}(20, 560:end)));
check('decoder: short output also aligned in the far range', c3 + 559 == 1200);

% decoder = direct convolution
f   = mc_filter(P.pulse(2), P);
z   = conv(double(v(5, :)), f.h.');
idx = (1:size(v, 2)) + f.shift;
ref = zeros(1, size(v, 2));
okI = idx <= numel(z);
ref(okI) = z(idx(okI));
check('decoder: FFT convolution = conv()', ...
      max(abs(double(each{2}(5, :)) - ref)) < 1e-4 * max(abs(ref)));

% coefficients given as a correlation reference
Pc = P;
Pc.decoder.coefForm = 'correlate';
for k = 1:2
    Pc.pulse(k).coef = conj(flipud(P.pulse(k).coef(:)));
end
yc = mc_decoder(v, Pc);
check('decoder: coefForm ''correlate'' = ''fir''', max(abs(yc(:) - y(:))) < 1e-4 * max(abs(y(:))));

% two frequencies + band filters, with the lag measured from the code or assumed centred
Pf = P0;
Pf.pulse(1).freqMHz = -0.5;
Pf.pulse(2).freqMHz = 0.5;
Pf.decoder.bandFilter = true;
[vf, ~, Pf] = mc_synthetic(Pf, S);
yf = abs(mc_decoder(vf, Pf));
[~, c1] = max(yf(20, 1:559));
[~, c2] = max(yf(20, 560:end));
check('decoder: separate frequencies + band filters', c1 == 300 && c2 + 559 == 1200);
Pn = Pf;
Pn.pulse(1).code = [];
Pn.pulse(2).code = [];
yn = abs(mc_decoder(vf, Pn));
[~, c1] = max(yn(20, 1:559));
[~, c2] = max(yn(20, 560:end));
check('decoder: centred-filter lag without the code', c1 == 300 && c2 + 559 == 1200);

% ---- canceler ------------------------------------------------------------------------------------
x = single(randn(50, 30) + 1j*randn(50, 30));
[yc, r] = mc_canceler(x, P0, mc_rows([], 50));
yr = filter([1 -2 1], 1, double(x), [], 1);
check('canceler = filter([1 -2 1])', max(max(abs(double(yc(3:end, :)) - yr(3:end, :)))) < 1e-4 ...
      && all(all(yc(1:2, :) == 0)) && isequal(find(~r.valid)', [1 2]));

% ---- Doppler FFT ----------------------------------------------------------------------------------
Pd = P0;
Pd.doppler.window = 'hamming';
Pd.doppler.hop = 5;
[D, r] = mc_doppler(x, Pd, mc_rows([], 50));
e = r.pulse(3);
w = mc_window('hamming', 16);
n = (0:15)';
F = zeros(16, 30);
for k = 1:16
    F(k, :) = sum(w .* double(x(e-15:e, :)) .* exp(-1j*2*pi*(k-1)*n/16), 1);
end
check('Doppler FFT = direct DFT', max(max(abs(squeeze(D(3, :, :)) - abs(F).'))) < 1e-3);
check('Doppler FFT rows end every hop pulses', isequal(r.pulse', 16:5:50));
t = exp(1j*2*pi*6/16*(0:49)') * ones(1, 3);                    % tone in bin 7
[D, ~] = mc_doppler(single(t), P0, mc_rows([], 50));
[~, b] = max(D(10, 1, :));
check('Doppler FFT: tone in its bin', b == 7);

% ---- non-coherent integration ----------------------------------------------------------------------
Sx = single(rand(30, 20, 4));
Pn = P0;
Pn.nci.length = 5;
[Y, r] = mc_nci(Sx, Pn, struct('pulse', (1:30)', 'seq', (1:30)', 'hasSeq', false, ...
                               'valid', true(30, 1)));
cs = cumsum(double(Sx), 1);
Yr = cs - [zeros(5, 20, 4); cs(1:end-5, :, :)];
check('NCI buffer = sum of the last L outputs', max(abs(double(Y(:)) - Yr(:))) < 1e-4 ...
      && isequal(find(~r.valid)', 1:4));

% ---- CFAR ------------------------------------------------------------------------------------------
A = single(rand(6, 80, 2)) + 0.1;
A(3, 40, 1) = 30;
A(5, 2, 2) = 30;
types = {'SO', 'CA', 'GO'};
edges = {'partial', 'skip'};
allOk = true;
for ti = 1:3
    for ei = 1:2
        Pq = P0;
        Pq.cfar.type = types{ti};
        Pq.cfar.edge = edges{ei};
        Pq.cfar.nRef = 6;
        Pq.cfar.nGuard = 2;
        Pq.cfar.thresholdDb = 6;
        [Dq, Tq] = mc_cfar(A, Pq);
        [Dr, Tr] = cfarLoop(A, Pq);
        allOk = allOk && isequal(Dq > 0, Dr > 0) && ...
                max(abs(Tq(isfinite(Tr)) - Tr(isfinite(Tr)))) < 1e-4;
    end
end
check('CFAR SO / CA / GO, partial / skip = direct loop', allOk);
Pq = P0;
Pq.cfar.nRef = 6;
Pq.cfar.nGuard = 2;
Dq = mc_cfar(A, Pq);
check('CFAR finds the two strong cells', Dq(3, 40, 1) == 30 && Dq(5, 2, 2) == 30);

% ---- comparison: known shift and gain are found ------------------------------------------------------
Pm = P0;
Pm.compare.plot = false;
y2 = mc_decoder(v, P);
rr = mc_rows(info, size(y2, 1));
refL = 2.5 * [y2(:, 8:end), zeros(size(y2, 1), 7, 'single')];         % ref(c) = 2.5 mine(c + 7)
fprintf('  (comparison output:)\n');
R = mc_compare(y2, rr, refL, info, 'test', Pm);
check('compare: shift +7 and gain +7.96 dB found', R.shift == 7 && abs(R.gainDb - 20*log10(2.5)) < 0.01 ...
      && R.coh > 0.9999);
R = mc_compare(abs(y2), rr, abs(y2).^2, info, 'test power', Pm);
check('compare: power against amplitude gives slope 2', abs(R.slope - 2) < 1e-3);

% ---- align check: the transmit delay of the short pulse is found -----------------------------------------
Pw = P;
Pw.pulse(1).txDelayUs = 238;                                        % 4 samples wrong
[~, eachW, filtW] = mc_decoder(v + single(0.01*(randn(size(v)) + 1j*randn(size(v)))), Pw);
A2 = mc_align_check(eachW, rr, y2, info, Pw, filtW);
check('align check: short pulse txDelayUs = 240 found', abs(A2(1).txDelayUs - 240) < 1e-9);

% ---- FFT rows aligned to the log rows --------------------------------------------------------------------
rc = mc_rows(info, size(v, 1));
logInfo.pulseSeq = info.pulseSeq(20:3:end);
e = mc_fft_ends(rc, logInfo, numel(logInfo.pulseSeq), P0);
check('FFT ends on the pulses of the log rows', isequal(e(:)', 20:3:size(v, 1)));

% ---- whole chain on the synthetic scene ------------------------------------------------------------------
[v, info, P, truth] = mc_synthetic(P0);
dec = mc_decoder(v, P);
[can, r] = mc_canceler(dec, P, mc_rows(info, size(v, 1)));
[dop, r] = mc_doppler(can, P, r);
[nci, r] = mc_nci(dop, P, r);
[M, B] = mc_max(mc_cfar(nci, P), P);
allOk = true;
for i = 1:numel(truth)
    hit = M(r.valid, truth(i).cell);
    bin = B(r.valid, truth(i).cell);
    allOk = allOk && all(hit > 0) && mode(double(bin)) == truth(i).bin;
end
check('chain: every target detected in every row, in its Doppler bin', allOk);

warning(quiet);
fprintf('%d passed, %d failed\n', nPass, nFail);
ok = nFail == 0;
end

function [D, T] = cfarLoop(Y, P)
% direct CFAR, one cell at a time
C = P.cfar;
alpha = 10^(C.thresholdDb / 20);
[nR, nC, nB] = size(Y);
D = zeros(size(Y));
T = nan(size(Y));
for b = 1:nB
    for i = 1:nR
        for c = 1:nC
            lead = (c - C.nGuard - C.nRef):(c - C.nGuard - 1);
            lag  = (c + C.nGuard + 1):(c + C.nGuard + C.nRef);
            lead = lead(lead >= 1);
            lag  = lag(lag <= nC);
            if strcmpi(C.edge, 'skip') && (numel(lead) < C.nRef || numel(lag) < C.nRef)
                continue
            end
            m1 = mean(double(Y(i, lead, b)));
            m2 = mean(double(Y(i, lag, b)));
            if isempty(lead), m1 = NaN; end
            if isempty(lag),  m2 = NaN; end
            switch C.type
                case 'SO', noise = min(m1, m2);
                case 'GO', noise = max(m1, m2);
                case 'CA', noise = mean(double(Y(i, [lead lag], b)));
            end
            T(i, c, b) = alpha * noise;
            if Y(i, c, b) > T(i, c, b)
                D(i, c, b) = Y(i, c, b);
            end
        end
    end
end
end
