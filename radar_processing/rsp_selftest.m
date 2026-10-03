function ok = rsp_selftest()
%RSP_SELFTEST Verify every block of the chain against reference computations.
%
%   ok = rsp_selftest()
%
%   1. matched filter : range alignment of both pulses and of the stitched output
%   2. canceler       : taps, clutter rejection, against a direct loop
%   3. Doppler FFT    : tone lands in the right bin, against a direct loop
%   4. integration    : against a direct loop
%   5. CFAR           : SO / CA / GO and every edge mode against a direct loop,
%                       threshold factor against the closed forms
%   6. full scenario  : targets found in the expected range cell and Doppler bin,
%                       false-alarm rate on noise close to the design Pfa
%   7. rsp_compare    : finds a known shift and gain

rng(7);
ok = true;
P  = rsp_default_params();
fprintf('rsp_selftest\n');

%% 1. Matched filter ----------------------------------------------------------------
M = 4;  R = 600;
cells = [30 300];                                  % blind zone and far range
S = struct('nPulses', M, 'nRange', R, 'noisePower', 1e-12, 'blankTx', false, 'seed', 1, ...
           'targets', struct('rangeCell', {cells(1), cells(2)}, 'fdNorm', {0, 0}, 'snrDb', {120, 120}), ...
           'clutter', []);
v = rsp_simulate(P, S);
[y, info] = rsp_matched_filter(v, P);
[~, pk1] = max(abs(info.each{info.shortPulse}(1, 1:info.switchCell)));
[~, pk2] = max(abs(info.each{info.longPulse}(1, info.switchCell+1:end)));
pk2 = pk2 + info.switchCell;
[~, pkS] = max(abs(y(1, 1:info.switchCell)));
ok = check(ok, 'matched filter: short pulse peak at its cell', pk1 == cells(1));
ok = check(ok, 'matched filter: long pulse peak at its cell',  pk2 == cells(2));
ok = check(ok, 'matched filter: stitched output keeps both',    pkS == cells(1) && ...
           abs(y(1, cells(2))) == abs(info.each{info.longPulse}(1, cells(2))));
% direct correlation of one row
h   = info.replica{info.longPulse};
d   = info.delay(info.longPulse);
row = double(v(2, :));
ref = zeros(1, R);
for r = 1:R
    n = r - 1 + d;
    m = 0:numel(h)-1;
    m = m(n + m <= R - 1);
    ref(r) = sum(row(n + m + 1) .* conj(h(m + 1)).');
end
ok = check(ok, 'matched filter: equals direct correlation', ...
           relErr(ref, info.each{info.longPulse}(2, :)) < 1e-5);
ok = check(ok, 'matched filter: noise gain is 1', abs(norm(h) - 1) < 1e-12);

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

%% 6. Full scenario --------------------------------------------------------------------
S = struct('nPulses', 350, 'nRange', 5469, 'noisePower', 1, 'blankTx', true, 'seed', 3, ...
           'targets', struct('rangeCell', {50, 1000, 2500, 3500}, ...
                             'fdNorm',    {0.25, 3/16, -0.31, 0.40}, ...
                             'snrDb',     {0, -10, -15, -12}), ...
           'clutter', struct('cnrDb', 40, 'cells', [1 1500], 'spreadNorm', 0.01));
v   = rsp_simulate(P, S);
out = rsp_chain(v, P);
L   = out.cfar.list;
for k = 1:numel(S.targets)
    t   = S.targets(k);
    bin = mod(round(t.fdNorm * P.fft.nfft), P.fft.nfft) + 1;
    hit = abs(L.range - t.rangeCell) <= 1 & abs(L.bin - bin) <= 1;
    ok  = check(ok, sprintf('scenario: target at cell %4d found in bin %2d (%d hits)', ...
                t.rangeCell, bin, nnz(hit)), nnz(hit) > 0);
end
% false alarms on noise only: beyond the clutter, away from the targets and
% from the far end (where the long pulse is only partly received)
lastFull = size(v, 2) - info.length(info.longPulse) - info.delay(info.longPulse);
inNoise  = L.range > 1600 & L.range < lastFull & ...
           all(abs(bsxfun(@minus, L.range, [S.targets.rangeCell])) > 30, 2);
nTest    = size(out.integral, 1) * size(out.integral, 3) * ...
           (lastFull - 1600 - 1 - 4*61);
pfa      = nnz(inNoise) / nTest;
ok = check(ok, sprintf('scenario: measured Pfa %.2g (design %.0g)', pfa, P.cfar.pfa), ...
           pfa < 5*P.cfar.pfa);

%% 7. rsp_compare ----------------------------------------------------------------------
ref  = out.decoder;
test = 0.5 * circshift(ref, [2 5]);
r = rsp_compare(ref, test, 'shifted decoder', 'quiet', true);
ok = check(ok, sprintf('compare: shift [%d %d], gain %.2f dB', r.lag, r.gainDb), ...
           isequal(r.lag, [2 5]) && abs(r.gainDb - 20*log10(2)) < 0.01 && r.corr > 0.999);
r = rsp_compare(out.cfar.map, out.cfar.det, 'cfar', 'quiet', true);
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
