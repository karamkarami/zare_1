function ok = rsp_selftest()
%RSP_SELFTEST Verify every block of the chain against reference computations.
%
%   ok = rsp_selftest()
%
%   1. codes / decoders : code properties, PSL, decoder forms, range alignment
%                         (matched, mismatched, LFM), against a direct loop
%   2. canceler         : taps, clutter rejection, against a direct loop
%   3. Doppler FFT      : tone lands in the right bin, against a direct loop
%   4. integration      : against a direct loop
%   5. CFAR             : SO / CA / GO and every edge mode against a direct loop,
%                         threshold factor against the closed forms
%   6. antenna, geometry, max stage after CFAR
%   7. block scenario   : targets (m, deg, m/s) in the expected cell and bin,
%                         false-alarm rate on noise close to the design Pfa
%   8. sector scan      : plots against truth, block size independence
%   9. rsp_compare      : finds a known shift and gain

rng(7);
ok = true;
P  = rsp_default_params();
fprintf('rsp_selftest\n');

%% 1. Codes, decoders and pulse compression -------------------------------------------
c = rsp_code('mls', 63);
a = real(ifft(abs(fft(c)).^2));
ok = check(ok, 'codes: MLS 63 periodic autocorrelation sidelobes = -1', all(abs(a(2:end) + 1) < 1e-9));
d = rsp_decoder(P.pulse(1), P.fs, 'noise');
ok = check(ok, sprintf('codes: Barker 13 matched PSL %.2f dB', d.pslDb), abs(d.pslDb + 22.28) < 0.01);
q = P.pulse(1);  q.decoder = rsp_code_mmf(q.code, 39);
d = rsp_decoder(q, P.fs, 'noise');
ok = check(ok, sprintf('codes: Barker 13 mismatched decoder PSL %.1f dB, loss %.2f dB', ...
           d.pslDb, -d.lossDb), d.pslDb < -35 && d.lossDb > -0.5);
q = P.pulse(1);  q.decoder = conj(fliplr(q.code));               % 'fir' matched taps
d1 = rsp_decoder(q, P.fs, 'noise');
q.decoder = q.code;  q.decoderForm = 'reference';               % same filter, other form
d2 = rsp_decoder(q, P.fs, 'noise');
d0 = rsp_decoder(P.pulse(1), P.fs, 'noise');
ok = check(ok, 'codes: decoder forms fir / reference / [] give the same filter', ...
           relErr(d0.h, d1.h) < 1e-12 && relErr(d0.h, d2.h) < 1e-12);

R = 2500;  cells = [200 1500];                                   % short and long pulse zones
configs = {'codes, matched', 'short pulse mismatched decoder', 'long pulse LFM + Taylor'};
for k = 1:numel(configs)
    Q = P;
    if k == 2
        Q.pulse(1).decoder = rsp_code_mmf(Q.pulse(1).code, 39);
    elseif k == 3
        Q.pulse(2).type = 'lfm';  Q.pulse(2).widthUs = 63;
        Q.pulse(2).window = 'taylor';  Q.pulse(2).windowParam = [4 -35];
    end
    Q.mf.keepEach = true;
    v = pointEchoes(Q, cells, 2, R);
    [y, info] = rsp_matched_filter(v, Q);
    [~, pk1] = max(abs(y(1, 1:info.switchCell)));
    [~, pk2] = max(abs(y(1, info.switchCell+1:end)));
    ok = check(ok, sprintf('decoder (%s): peaks at cells %d and %d', configs{k}, ...
               pk1, pk2 + info.switchCell), pk1 == cells(1) && pk2 + info.switchCell == cells(2));
end
% direct correlation of one row (long pulse, with alignment lag)
h   = info.replica{info.longPulse};
n0  = info.delay(info.longPulse) + info.dec{info.longPulse}.lag;
row = double(v(2, :));
ref = zeros(1, R);
for r = 1:R
    for m = 0:numel(h)-1
        n = r - 1 + n0 + m;
        if n >= 0 && n <= R - 1
            ref(r) = ref(r) + row(n + 1) * conj(h(m + 1));
        end
    end
end
ok = check(ok, 'decoder: FFT correlation equals direct loop', ...
           relErr(ref, info.each{info.longPulse}(2, :)) < 1e-5);
ok = check(ok, 'decoder: noise gain is 1', abs(norm(h) - 1) < 1e-12);

%% 2. Canceler -------------------------------------------------------------------------
X = complex(randn(40, 30), randn(40, 30));
[Y, idx, w] = rsp_canceler(X, P.canceler);
Yd = zeros(size(X));
for m = 3:size(X, 1)
    Yd(m, :) = X(m, :) - 2*X(m-1, :) + X(m-2, :);
end
ok = check(ok, 'canceler: 3-pulse taps [1 -2 1]', isequal(w, [1 -2 1]));
ok = check(ok, 'canceler: equals direct loop', relErr(Yd, Y) < 1e-12 && numel(idx) == 40);
Yc = rsp_canceler(ones(10, 5), P.canceler);
ok = check(ok, 'canceler: rejects a constant (zero Doppler)', max(abs(Yc(:))) < 1e-12);
C2 = P.canceler;  C2.output = 'valid';
[Yv, idxv] = rsp_canceler(X, C2);
ok = check(ok, 'canceler: valid output drops 2 rows', size(Yv, 1) == 38 && idxv(1) == 3);

%% 3. Doppler FFT ----------------------------------------------------------------------
F = P.fft;  F.hop = 3;
bin = 6;
X = repmat(exp(1j*2*pi*(bin-1)/F.nfft*(0:49)'), 1, 7);
[D, idx] = rsp_doppler_fft(X, F);
[~, b] = max(abs(squeeze(D(2, 4, :))));
ok = check(ok, 'Doppler FFT: tone in the expected bin', b == bin);
wF = rsp_window(F.window, F.nPulses);
f  = 3;  s = (f-1)*F.hop + (1:F.nPulses);
Dd = fft(X(s, :) .* wF, F.nfft, 1).';
ok = check(ok, 'Doppler FFT: equals direct FFT', relErr(Dd, squeeze(D(f, :, :))) < 1e-12 && ...
           idx(f) == s(end));
F.chunk = 2;
ok = check(ok, 'Doppler FFT: independent of chunk size', relErr(D, rsp_doppler_fft(X, F)) < 1e-14);

%% 4. Integration ----------------------------------------------------------------------
D = complex(randn(20, 9, 4), randn(20, 9, 4));
C = P.nci;
[S4, idx] = rsp_nci(D, C);
n  = C.nFrames;
Sd = zeros(20 - n + 1, 9, 4);
for f = n:20
    Sd(f - n + 1, :, :) = mean(abs(D(f-n+1:f, :, :)).^2, 1);
end
ok = check(ok, 'integration: equals direct mean of |x|^2', relErr(Sd, S4) < 1e-12 && idx(1) == n);

%% 5. CFAR -----------------------------------------------------------------------------
X = -log(rand(6, 120, 3));                           % exponential noise
X(:, 60, 2) = 200;                                   % one target
for type = {'SO', 'CA', 'GO'}
    for edge = {'oneSided', 'partial', 'none'}
        C = P.cfar;
        C.type = type{1};  C.edge = edge{1};  C.nRef = 8;  C.nGuard = 2;
        C.thresholdMode = 'factor';  C.factorDb = 10;
        c = rsp_cfar(X, C, 'square');
        T = directCfar(X, C, 10);
        ok = check(ok, sprintf('CFAR %s / %-8s: equals direct loop', type{1}, edge{1}), ...
                   isequal(c.det, X > T) && relErr(T(isfinite(T)), double(c.threshold(isfinite(T)))) < 1e-6);
    end
end
ok = check(ok, 'CFAR: target detected', c.det(1, 60, 2));
aCA = rsp_cfar_factor('CA', 16, 1, 1e-6);
aSO = rsp_cfar_factor('SO', 5, 1, 1e-6);
T   = aSO/5;  kk = 0:4;
pSO = 2*(2+T)^(-5)*sum(arrayfun(@(k) nchoosek(4+k, k), kk).*(2+T).^(-kk));
ok = check(ok, 'CFAR factor: CA closed form', abs(aCA - 32*(1e-6^(-1/32) - 1)) < 1e-6*aCA);
ok = check(ok, 'CFAR factor: SO closed form (Weiss)', abs(pSO - 1e-6) < 1e-9);

%% 6. Antenna, geometry, output stage --------------------------------------------------
G = 10*log10(rsp_antenna_pattern([0 2 -2 6.4 180 362], P.antenna));
ok = check(ok, 'antenna: 3 dB at +-2 deg, sidelobe -25 dB, back lobe <= -40 dB, wraps', ...
           abs(G(2) + 3.01) < 0.01 && abs(G(3) + 3.01) < 0.01 && abs(G(4) + 25) < 0.01 && ...
           G(5) <= -40 && abs(G(6) - G(2)) < 1e-9);
g = rsp_geometry(P, 10, (1:3)');
ok = check(ok, sprintf('geometry: %.1f m cells, %.3f deg/pulse, +-%.1f m/s unambiguous', ...
           g.cellM, g.degPerPulse, g.vUnambMps/2), abs(g.cellM - 24.98) < 0.01 && ...
           abs(g.degPerPulse - 0.036) < 1e-12 && abs(g.binVelMps(2) - g.vUnambMps/16) < 1e-9);
map = single(rand(5, 7, 16) .* (rand(5, 7, 16) > 0.7));
[v, b] = rsp_cfar_max(map, [2:8 10]);
[vd, kd] = max(map(:, :, [2:8 10]), [], 3);
bd = [2:8 10];  bd = uint8(bd(kd));  bd(vd <= 0) = 0;
ok = check(ok, 'output: max over selected bins after CFAR', isequal(v, vd) && isequal(b, bd));
X = -log(rand(4, 200, 3));  X(:, 100, :) = 300;
C = P.cfar;  C.thresholdMode = 'factor';  C.nRef = 16;  C.nGuard = 2;  C.mapValue = 'snr';
c = rsp_cfar(X, C, 'square');
nz = c.threshold(c.det) / 10^(C.factorDb/10);
ok = check(ok, 'CFAR: snr output = value / noise estimate', ...
           relErr(double(X(c.det)) ./ double(nz), c.map(c.det)) < 1e-5);

%% 7. Block scenario (350 pulses, targets in m / deg / m/s) ---------------------------------
S = struct('nRange', 5469, 'noisePower', 1, 'blankTx', true, 'seed', 3);
S.targets = struct('rangeM', {6000, 40000, 90000}, 'azDeg', {5, 6, 7}, ...
                   'velocityMps', {40, -75, 200}, 'snrDb', {-10, -15, -18});
S.clutter = struct('cnrDb', 20, 'maxRangeM', 25000, 'sigmaVMps', 0.5, 'textureDb', 3);
v   = rsp_simulate(P, S, (1:350)');
out = rsp_chain(v, P);
tr  = rsp_truth(P, S);
for k = 1:numel(tr)
    [f, c] = find(out.max.value > 0);
    b   = out.max.bin(sub2ind(size(out.max.value), f, c));
    hit = abs(c - tr(k).cell) <= 1 & abs(double(b) - tr(k).bin) <= 1;
    ok  = check(ok, sprintf('block: target %5.0f m, %4.0f m/s -> cell %d, bin %d (%d hits)', ...
                tr(k).rangeM, tr(k).velocityMps, tr(k).cell, tr(k).bin, nnz(hit)), nnz(hit) > 0);
end
% false alarms on noise only
S0 = S;  S0.targets = S.targets([]);  S0.clutter = [];  S0.seed = 4;
o0 = rsp_chain(rsp_simulate(P, S0, (1:350)'), P);
vc = o0.mf.validCells;
pfa = nnz(o0.cfar.det) / (size(o0.integral, 1) * size(o0.integral, 3) * (vc(2) - vc(1) + 1));
ok = check(ok, sprintf('block: measured Pfa on noise %.2g (design %.0g)', pfa, P.cfar.pfa), ...
           pfa < 5*P.cfar.pfa);

%% 8. Sector scan and plot extraction ---------------------------------------------------------
S.targets = struct('rangeM', {30000, 75000}, 'azDeg', {8, 14}, ...
                   'velocityMps', {-60, 120}, 'snrDb', {-15, -18});
scan = rsp_scan(P, S, 'degrees', 22, 'display', false, 'quiet', true);
rep  = scan.report;
ok = check(ok, sprintf('scan: %d plots for %d targets', numel(scan.plots.azDeg), numel(rep)), ...
           numel(scan.plots.azDeg) == numel(rep) && all([rep.found]));
if all([rep.found])
    ok = check(ok, sprintf('scan: errors az %.2f deg, range %.0f m, velocity %.1f m/s', ...
               max(abs([rep.dAzDeg])), max(abs([rep.dRangeM])), max(abs([rep.dVelMps]))), ...
               max(abs([rep.dAzDeg])) < 0.3 && max(abs([rep.dRangeM])) < 30 && ...
               max(abs([rep.dVelMps])) <= g.vUnambMps/16);           % one Doppler bin
end
% block processing equals one-shot processing
S.targets = S.targets(1);  S.clutter = [];  S.seed = 9;
Q = P;  Q.ppi.blockPulses = 100;
s1 = rsp_scan(Q, S, 'degrees', 12, 'display', false, 'quiet', true);
Q.ppi.blockPulses = 1000;
s2 = rsp_scan(Q, S, 'degrees', 12, 'display', false, 'quiet', true);
[d1, i1] = sortrows([s1.det.pulse s1.det.cell]);
[d2, i2] = sortrows([s2.det.pulse s2.det.cell]);
ok = check(ok, sprintf('scan: block size does not change the output (%d detections)', ...
           size(d1, 1)), isequal(d1, d2) && relErr(s1.det.value(i1), s2.det.value(i2)) < 1e-5);

%% 9. rsp_compare ----------------------------------------------------------------------
ref  = out.decoder .* single(1 + 3*rand(size(out.decoder, 1), 1));   % structure along pulses
test = 0.5 * circshift(ref, [2 5]);
r = rsp_compare(ref, test, 'shifted decoder', 'quiet', true);
ok = check(ok, sprintf('compare: shift [%d %d], gain %.2f dB', r.lag, r.gainDb), ...
           isequal(r.lag, [2 5]) && abs(r.gainDb - 20*log10(2)) < 0.01 && r.corr > 0.999);
r = rsp_compare(out.cfar.map, out.cfar.det, 'cfar', 'quiet', true, 'align', false);
ok = check(ok, 'compare: identical detection maps', r.missed == 0 && r.extra == 0);

if ok
    fprintf('ALL TESTS PASSED\n');
else
    fprintf('SOME TESTS FAILED\n');
end
end

% ======================================================================================
function ok = check(ok, label, pass)
if pass
    fprintf('  [ ok ] %s\n', label);
else
    fprintf('  [FAIL] %s\n', label);
end
ok = ok && pass;
end

function v = pointEchoes(P, cells, M, R)
% Noise-free video: a point echo of the whole transmission at each cell.
v = zeros(M, R);
for k = 1:numel(P.pulse)
    s = rsp_waveform(P.pulse(k), P.fs);
    d = round(P.pulse(k).delayUs * P.fs);
    for c = cells
        n = (c - 1) + d + (0:numel(s)-1);
        in = n < R;
        v(:, n(in) + 1) = v(:, n(in) + 1) + repmat(s(in).', M, 1);
    end
end
v = single(v);
end

function e = relErr(a, b)
a = double(a(:));
b = double(b(:));
e = norm(a - b) / max(norm(a), realmin);
end

function T = directCfar(X, C, factorDb)
% Cell-by-cell CFAR threshold (reference for rsp_cfar).
[nF, R, nB] = size(X);
N = C.nRef;  G = C.nGuard;  a = 10^(factorDb/10);
T = inf(nF, R, nB);
for r = 1:R
    lead = max(r-G-N, 1):r-G-1;
    lag  = r+G+1:min(r+G+N, R);
    fl = numel(lead) == N;  fg = numel(lag) == N;
    switch lower(C.edge)
        case 'onesided'
            uL = fl || (~fg && ~isempty(lead));
            uG = fg || (~fl && ~isempty(lag));
        case 'partial'
            uL = ~isempty(lead);  uG = ~isempty(lag);
        case 'none'
            uL = fl && fg;  uG = uL;
    end
    if ~(uL || uG)
        continue
    end
    for f = 1:nF
        for b = 1:nB
            mL = NaN;  mG = NaN;
            if uL, mL = mean(X(f, lead, b)); end
            if uG, mG = mean(X(f, lag, b));  end
            switch C.type
                case 'SO', z = min(mL, mG);
                case 'GO', z = max(mL, mG);
                case 'CA', z = [mL mG];  z = mean(z(~isnan(z)));
            end
            T(f, r, b) = a * z;
        end
    end
end
end
