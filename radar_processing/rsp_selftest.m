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
%   7. simulator        : echo against an independent computation, Doppler
%                         phase, noise and clutter power levels
%   8. scenario         : the 3 targets (10 km 20 m/s 40 deg, 50 km 60 m/s 20 deg,
%                         80 km 100 m/s 80 deg) found with the right range,
%                         azimuth and velocity, CFAR margin as budgeted; Pfa
%                         on noise; block size independence
%   9. rsp_compare      : finds a known shift and gain
%  10. antenna patterns : sinc, a pattern loaded from a .csv / .mat file
%                         (beamwidth measured, values interpolated)
%  11. terrain clutter  : power and spectral spread of a forest layer, mountain
%                         shadow on clutter and on targets, region map
%  12. range ambiguity  : echo beyond c/2PRF at its folded range with the true phase
%  13. canceler / FFT / carrier / PRF : 2- and 4-pulse taps, 32-pulse FFT; a target
%                         through the whole chain with a 2-pulse canceler, with a
%                         32-pulse FFT and at PRF 1500 / 1300 MHz (velocity folding)
%                         found where expected, margins as budgeted
%  14. scan schedule    : timeline continuity (time, azimuth, gap), parameter
%                         paths, false alarm accounting per region
%  15. saving           : results.mat and CSV files read back unchanged

rng(7);
ok = true;
P  = rsp_default_params();
fprintf('rsp_selftest\n');

%% 1. Codes, decoders and pulse compression -------------------------------------------
c = rsp_code('mls', 63);
a = real(ifft(abs(fft(c)).^2));
ok = check(ok, 'codes: MLS 63 periodic autocorrelation sidelobes = -1', all(abs(a(2:end) + 1) < 1e-9));
q0 = P.pulse(1);  q0.decoder = [];                               % matched filter
d = rsp_decoder(q0, P.fs, 'noise');
ok = check(ok, sprintf('codes: Barker 13 matched PSL %.2f dB', d.pslDb), abs(d.pslDb + 22.28) < 0.01);
q = P.pulse(1);  q.decoder = rsp_code_mmf(q.code, 39);
d = rsp_decoder(q, P.fs, 'noise');
ok = check(ok, sprintf('codes: Barker 13 mismatched decoder PSL %.1f dB, loss %.2f dB', ...
           d.pslDb, -d.lossDb), d.pslDb < -35 && d.lossDb > -0.5);
q = q0;  q.decoder = conj(fliplr(q.code));                       % 'fir' matched taps
d1 = rsp_decoder(q, P.fs, 'noise');
q.decoder = q.code;  q.decoderForm = 'reference';               % same filter, other form
d2 = rsp_decoder(q, P.fs, 'noise');
d0 = rsp_decoder(q0, P.fs, 'noise');
ok = check(ok, 'codes: decoder forms fir / reference / [] give the same filter', ...
           relErr(d0.h, d1.h) < 1e-12 && relErr(d0.h, d2.h) < 1e-12);

R = 2500;  cells = [200 1500];                                   % short and long pulse zones
configs = {'Barker mismatched + LFM, channel filters', 'Barker matched, no channel filter', ...
           'long pulse MLS 63 code'};
for k = 1:numel(configs)
    Q = P;
    if k == 2
        Q.pulse(1).decoder = [];  Q.mf.channelFilter = false;
    elseif k == 3
        Q.pulse(2).type = 'code';  Q.pulse(2).code = rsp_code('mls', 63);
        Q.pulse(2).window = 'none';  Q.pulse(1).delayUs = 63;
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
C = P.nci;  C.nFrames = 4;
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
mx = rsp_cfar_max(map, [2:8 10]);  v = mx.value;  b = mx.bin;
[vd, kd] = max(map(:, :, [2:8 10]), [], 3);
bd = [2:8 10];  bd = uint8(bd(kd));  bd(vd <= 0) = 0;
ok = check(ok, 'output: max over selected bins after CFAR', isequal(v, vd) && isequal(b, bd));
X = -log(rand(4, 200, 3));  X(:, 100, :) = 300;
C = P.cfar;  C.thresholdMode = 'factor';  C.nRef = 16;  C.nGuard = 2;  C.mapValue = 'snr';
c = rsp_cfar(X, C, 'square');
nz = c.threshold(c.det) / 10^(C.factorDb/10);
ok = check(ok, 'CFAR: snr output = value / noise estimate', ...
           relErr(double(X(c.det)) ./ double(nz), c.map(c.det)) < 1e-5);

%% 7. Simulator against an independent computation --------------------------------------
c0 = 299792458;
Q  = P;  Q.antenna.beamwidthDeg = 4;
g  = rsp_geometry(Q, 4000);
cell0 = 1500;                                       % integer delay: exact comparison
Rt = (cell0 - 1) * g.cellM;
St = struct('nRange', 4000, 'noiseDb', -300, 'blankTx', false, 'seed', 1, 'clutter', []);
St.targets = struct('rangeM', Rt, 'azDeg', 0, 'velocityMps', 0, 'powerDb', 10);
v  = double(rsp_simulate(Q, St, 1));                % pulse 1: antenna on the target
tx = pointEchoes(Q, 1, 1, 4000);                    % whole transmission at cell 1
ref = circshift(double(tx), [0 cell0 - 1]) * sqrt(10) * exp(-1j*2*pi*Q.radar.fcMHz*2*Rt/c0*1e6);
ok = check(ok, 'simulator: echo = delayed transmission x carrier phase x amplitude', ...
           relErr(ref, v) < 1e-6);                     % single precision video
St.targets.velocityMps = 30;                        % Doppler and its sign
v  = rsp_simulate(Q, St, (1:8)');
[~, n] = max(abs(v(1, :)));
dph = angle(sum(v(2:end, n) .* conj(v(1:end-1, n)))) / (2*pi);
ok = check(ok, sprintf('simulator: Doppler %.4f cycles/pulse (expected 2v/lambda/PRF = %.4f)', ...
           dph, 2*30/g.lambdaM/Q.radar.prfHz), abs(dph - 2*30/g.lambdaM/Q.radar.prfHz) < 1e-3);
St.targets.velocityMps = 0;  St.noiseDb = 0;  St.targets.powerDb = -300;    % levels
St.clutter = struct('powerDb', -50, 'maxRangeM', 60000, 'sigmaVMps', 0.5, 'textureDb', 0);
v  = double(rsp_simulate(Q, St, (1:200)'));
x  = v(:, 3500:4000);                               % beyond the clutter
pn = 10*log10(mean(abs(x(:)).^2));
St.noiseDb = -300;
v  = double(rsp_simulate(Q, St, (1:200)'));
x  = v(:, 900:1500);                                % full transmission overlap, before the fade
pc = 10*log10(mean(abs(x(:)).^2));
ok = check(ok, sprintf('simulator: noise %.2f dB (0), clutter %.2f dB (-50)', pn, pc), ...
           abs(pn) < 0.1 && abs(pc + 50) < 0.5);

%% 8. Your scenario: 3 targets through the whole chain ---------------------------------------
S = struct('nRange', 5469, 'noiseDb', 0, 'blankTx', true, 'seed', 3);
S.targets = struct('rangeM', {10000, 50000, 80000}, 'azDeg', {40, 20, 80}, ...
                   'velocityMps', {20, 60, 100}, 'powerDb', {10, 5, 0});
S.clutter = struct('powerDb', -50, 'maxRangeM', 25000, 'sigmaVMps', 0.5, 'textureDb', 3);
scan = rsp_scan(P, S, 'degrees', 95, 'display', false, 'quiet', true);
rep  = scan.report;
bud  = scan.budget;
g    = rsp_geometry(P, S.nRange);
ok = check(ok, sprintf('scenario: %d plots for %d targets, all found', numel(scan.plots.azDeg), ...
           numel(rep)), numel(scan.plots.azDeg) == numel(rep) && all([rep.found]));
if all([rep.found])
    ok = check(ok, sprintf('scenario: errors az %.2f deg, range %.0f m, velocity %.1f m/s', ...
               max(abs([rep.dAzDeg])), max(abs([rep.dRangeM])), max(abs([rep.dVelMps]))), ...
               max(abs([rep.dAzDeg])) < 0.2 && max(abs([rep.dRangeM])) < 25 && ...
               max(abs([rep.dVelMps])) < 3);
    for k = 1:numel(rep)
        ok = check(ok, sprintf(['scenario: target %d passes the CFAR by %.1f dB ' ...
                   '(budget %.1f dB)'], rep(k).target, rep(k).marginDb, ...
                   bud(rep(k).target).marginDb), rep(k).marginDb > 10 && ...
                   abs(rep(k).marginDb - bud(rep(k).target).marginDb) < 4);
    end
end
% false alarms on noise only
S0 = S;  S0.targets = S.targets([]);  S0.clutter = [];  S0.seed = 4;
o0 = rsp_chain(rsp_simulate(P, S0, (1:400)'), P);
vc = o0.mf.validCells;
pfa = nnz(o0.cfar.det) / (size(o0.integral, 1) * size(o0.integral, 3) * (vc(2) - vc(1) + 1));
ok = check(ok, sprintf('noise only: measured Pfa %.2g (design %.0g)', pfa, P.cfar.pfa), ...
           pfa < 5*P.cfar.pfa);
out = o0;
% block processing equals one-shot processing
S.targets = S.targets(2);  S.clutter = [];  S.seed = 9;
S.targets.azDeg = 6;
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

%% 10. Antenna patterns ---------------------------------------------------------------
A = P.antenna;  A.type = 'sinc';
G = 10*log10(rsp_antenna_pattern([0 A.beamwidthDeg/2 -A.beamwidthDeg/2], A));
th = (1.2:0.001:2.2) * A.beamwidthDeg;                            % first sidelobe region
sl = max(10*log10(rsp_antenna_pattern(th, A)));
ok = check(ok, sprintf('antenna sinc: 3 dB at +-bw/2, first sidelobe %.2f dB (-13.26)', sl), ...
           abs(G(2) + 3.01) < 0.02 && abs(G(3) + 3.01) < 0.02 && abs(sl + 13.26) < 0.02);
here = fileparts(mfilename('fullpath'));
A  = rsp_antenna_load(P.antenna, fullfile(here, 'patterns', 'example_pattern.csv'));
raw = dlmread(fullfile(here, 'patterns', 'example_pattern.csv'), ',', 1, 0);
Gt = 10*log10(rsp_antenna_pattern(raw(:, 1), A));
G3 = 10*log10(rsp_antenna_pattern(A.beamwidthDeg/2 * [-1 1], A));
ok = check(ok, sprintf(['antenna file: %d points, values reproduced, measured beamwidth %.2f deg ' ...
           '(3 dB edges %.2f / %.2f dB)'], size(raw, 1), A.beamwidthDeg, G3), ...
           strcmp(A.type, 'table') && max(abs(Gt - (raw(:, 2) - max(raw(:, 2))))) < 1e-9 && ...
           abs(A.beamwidthDeg - 4.02) < 0.05 && all(abs(G3 + 3) < 0.05));
tmp = [tempname() '.mat'];
angleDeg = 0:0.5:180;  gainDb = -12 * (angleDeg / 4).^2;  gainDb = max(gainDb, -45); %#ok<NASGU>
save(tmp, 'angleDeg', 'gainDb');
A = rsp_antenna_load(P.antenna, tmp);
delete(tmp);
ok = check(ok, sprintf('antenna .mat one-sided table: symmetric, beamwidth %.3f deg (4)', ...
           A.beamwidthDeg), abs(A.beamwidthDeg - 4) < 0.01 && ...
           abs(rsp_antenna_pattern(-3, A) - rsp_antenna_pattern(3, A)) < 1e-12);

%% 11. Terrain clutter --------------------------------------------------------------------
g  = rsp_geometry(P, 10);
St = struct('nRange', 1200, 'noiseDb', -300, 'blankTx', false, 'seed', 5);
St.targets = struct('rangeM', {}, 'azDeg', {}, 'velocityMps', {}, 'powerDb', {});
St.clutter = struct('type', 'forest', 'azDeg', [0 360], 'rangeM', [0 60000], 'powerDb', 20, ...
                    'textureDb', 0, 'edgeM', 0, 'dcFraction', 0.5, 'sigmaVMps', 0.6);
v  = double(rsp_simulate(P, St, (1:2000)'));
x  = v(:, 700:1150);
pc = 10*log10(mean(abs(x(:)).^2));
a  = x(2:end, :);  b = x(1:end-1, :);
r1 = abs(sum(a(:) .* conj(b(:)))) / sqrt(sum(abs(a(:)).^2) * sum(abs(b(:)).^2));
sv = sqrt(2*(1 - r1)) / (2*pi) * g.lambdaM * P.radar.prfHz / 2;   % spectral std [m/s]
vs = 0.265 * g.degPerPulse / P.antenna.beamwidthDeg * g.lambdaM * P.radar.prfHz / 2;
ve = sqrt(0.5 * 0.6^2 + vs^2);
ok = check(ok, sprintf(['clutter forest: power %.2f dB (20), spectral spread %.3f m/s ' ...
           '(%.3f: half static, half 0.6 m/s wind, + %.2f m/s scanning)'], pc, sv, ve, vs), ...
           abs(pc - 20) < 0.5 && abs(sv - ve) < 0.1 * ve);
% mountain shadow on the ground clutter behind it
St.clutter = struct('type', {'ground', 'mountain'}, 'azDeg', {[0 360], [100 140]}, ...
                    'rangeM', {[0 60000], [15000 20000]}, 'powerDb', {20, 20}, ...
                    'textureDb', {0, 0}, 'edgeM', {0, 0});
St.nRange = 2000;
pw = zeros(1, 2);
Q  = P;
for i = 1:2
    Q.radar.azStartDeg = 120 + 180 * (i - 1);                  % mountain / open side
    v = double(rsp_simulate(Q, St, (1:100)'));
    x = v(:, 1600:1990);                                        % 40 .. 50 km: the video spans
                                                                % the 17 km transmission
    pw(i) = 10*log10(mean(abs(x(:)).^2));
end
ok = check(ok, sprintf('clutter: shadow behind the mountain %.1f dB (open side %.1f dB)', ...
           pw(1) - pw(2), pw(2)), pw(1) - pw(2) < -35 && abs(pw(2) - 20) < 0.5);
% shadow on a target, truth and region map
St.clutter(2).powerDb = -300;                                    % shadow only
St.targets = struct('rangeM', 30000, 'azDeg', 120, 'velocityMps', 0, 'powerDb', 0);
St.noiseDb = -300;  St.clutter(1).powerDb = -300;
Q.radar.azStartDeg = 120;
v1 = max(abs(double(rsp_simulate(Q, St, 1))));
tr = rsp_truth(Q, St);
St.clutter = [];
v0 = max(abs(double(rsp_simulate(Q, St, 1))));
L  = rsp_clutter_types(struct('type', {'ground', 'mountain', 'forest'}));
[id, nm] = rsp_clutter_region(L, [10 120 120 225 120], [5000 47000 70000 20000 30000]);
ok = check(ok, sprintf(['clutter: target behind the mountain %.1f dB (truth %.1f dB); ' ...
           'regions %s / %s / %s / %s / %s'], 20*log10(v1/v0), tr.shadowDb, nm{id + 1}), ...
           abs(20*log10(v1/v0) + 40) < 0.1 && abs(tr.shadowDb + 40) < 0.01 && ...
           isequal(id, [1 2 4 3 0]));

%% 12. Range ambiguity: echo of an earlier pulse ---------------------------------------
Q  = P;  Q.radar.prfHz = 2000;  Q.antenna.type = 'model';
g  = rsp_geometry(Q, 10);
cell0 = 1500;
Rt = g.rangeUnambM + (cell0 - 1) * g.cellM;                     % beyond c/2PRF = 75 km
St = struct('nRange', 2900, 'noiseDb', -300, 'blankTx', false, 'seed', 1, 'clutter', []);
St.targets = struct('rangeM', Rt, 'azDeg', 0, 'velocityMps', 0, 'powerDb', 10);
v  = double(rsp_simulate(Q, St, 1));
tx = pointEchoes(Q, 1, 1, 2900);
ref = circshift(double(tx), [0 cell0 - 1]) * sqrt(10) * exp(-1j*2*pi*Q.radar.fcMHz*2*Rt/c0*1e6);
tr = rsp_truth(Q, St);
ok = check(ok, sprintf(['range ambiguity: target at %.1f km (beyond %.1f km) appears at ' ...
           '%.1f km with the true carrier phase'], Rt/1e3, g.rangeUnambM/1e3, ...
           tr.apparentRangeM/1e3), relErr(ref, v) < 1e-6 && tr.ambiguous && tr.cell == cell0);

%% 13. Canceler order, FFT length, carrier and PRF through the whole chain ------------------
C = P.canceler;  C.order = 2;
X = complex(randn(30, 8), randn(30, 8));
[Y, ~, w] = rsp_canceler(X, C);
ok = check(ok, 'canceler: 2-pulse taps [1 -1], equals x(m) - x(m-1)', isequal(w, [1 -1]) && ...
           relErr([zeros(1, 8); X(2:end, :) - X(1:end-1, :)], Y) < 1e-12);
C.order = 4;
ok = check(ok, 'canceler: 4-pulse taps [1 -3 3 -1]', isequal(rsp_canceler_taps(C), [1 -3 3 -1]));
F = P.fft;  F.nPulses = 32;  F.nfft = 32;
X = repmat(exp(1j*2*pi*(10-1)/32*(0:79)'), 1, 3);
[D, idx] = rsp_doppler_fft(X, F);
[~, b] = max(abs(squeeze(D(5, 2, :))));
s5 = (5-1)*F.hop + (1:32);
ok = check(ok, 'Doppler FFT 32 pulses: tone in bin 10, equals direct FFT', b == 10 && ...
           relErr(fft(X(s5, :) .* rsp_window(F.window, 32), 32, 1).', squeeze(D(5, :, :))) < 1e-12 ...
           && idx(5) == s5(end));
St = struct('nRange', 3000, 'noiseDb', 0, 'blankTx', true, 'seed', 11, 'clutter', []);
St.targets = struct('rangeM', 50000, 'azDeg', 6, 'velocityMps', 100, 'powerDb', 0);
cfgs = {'2-pulse canceler', {'canceler.order', 2}; ...
        'FFT 32 pulses, NCI 16', {'fft.nPulses', 32, 'fft.nfft', 32, 'nci.nFrames', 16}; ...
        'PRF 1500 Hz, 1300 MHz', {'radar.prfHz', 1500, 'radar.fcMHz', 1300}};
for k = 1:size(cfgs, 1)
    sc = struct('name', cfgs{k, 1}, 'degrees', 12, 'set', cfgs(k, 2));
    rn = rsp_schedule(P, St, sc, 'display', false, 'quiet', true);
    s1 = rn.scans(1);
    r  = s1.report;
    gk = rsp_geometry(s1.P, 10);
    if r.found
        ok = check(ok, sprintf(['%s: target found, errors %.2f deg %.0f m %.1f m/s ' ...
                   '(measures %.1f of 100 m/s, +-%.0f), margin %.1f dB (budget %.1f)'], ...
                   cfgs{k, 1}, r.dAzDeg, r.dRangeM, r.dVelMps, r.plotVelocityMps, ...
                   gk.vUnambMps/2, r.marginDb, s1.budget.marginDb), abs(r.dAzDeg) < 0.2 && ...
                   abs(r.dRangeM) < 25 && abs(r.dVelMps) < 3 && ...
                   abs(r.marginDb - s1.budget.marginDb) < 4);
    else
        ok = check(ok, sprintf('%s: target found', cfgs{k, 1}), false);
    end
end
ok = check(ok, sprintf('carrier / PRF: 100 m/s at 1300 MHz, PRF 1500 folds to %.1f m/s', ...
           s1.truth.foldedMps), abs(s1.truth.foldedMps - (100 - 2*86.48)) < 0.1);

%% 14. Scan schedule -----------------------------------------------------------------------
St = struct('nRange', 3000, 'noiseDb', 0, 'blankTx', true, 'seed', 21);
St.targets = struct('rangeM', {}, 'azDeg', {}, 'velocityMps', {}, 'powerDb', {});
St.clutter = struct('type', 'forest', 'azDeg', [0 30], 'rangeM', [10000 40000], 'powerDb', 20);
sc = struct('name', {'a', 'b'}, 'degrees', {5, 5}, 'gapS', {0, 10}, ...
            'set', {{'radar.prfHz', 1200}, {'pulse(1).fcMHz', -1, 'canceler.order', 2}});
rn = rsp_schedule(P, St, sc, 'display', false, 'quiet', true);
tl = rn.timeline;
s2 = rn.scans(2);
ok = check(ok, sprintf(['schedule: entry 2 starts at %.3f s = end %.3f s + gap 10 s, ' ...
           'az %.2f deg = %.2f + one idle turn'], tl.startS(2), tl.endS(1), tl.startAzDeg(2), ...
           tl.endAzDeg(1)), abs(tl.startS(2) - tl.endS(1) - 10) < 1e-9 && ...
           abs(mod(tl.startAzDeg(2) - tl.endAzDeg(1) + 180, 360) - 180) < 1e-9 && ...
           tl.prfHz(1) == 1200 && tl.prfHz(2) == P.radar.prfHz && ...
           s2.P.pulse(1).fcMHz == -1 && tl.canceler(2) == 2 && tl.canceler(1) == 3 && ...
           abs(s2.P.radar.startTimeS - tl.startS(2)) < 1e-12);
fa = s2.falseAlarms;
nTested = numel(s2.azDeg) * (diff(s2.validCells) + 1);
ok = check(ok, sprintf(['false alarms: %d cells counted in regions %s (= %d tested), ' ...
           'totals consistent'], fa.cells(end), strjoin(fa.names(fa.cells > 0), '/'), nTested), ...
           fa.cells(end) == nTested && fa.falseDet(end) == nnz(fa.isFalse) && ...
           fa.falseDet(end) == sum(fa.falseDet(1:end-1)) && fa.cells(2) > 0);

%% 15. Saving ------------------------------------------------------------------------------
dirOut = tempname();
files = evalc('rsp_save(rn, dirOut, ''figures'', false);');     %#ok<NASGU>
sv = load(fullfile(dirOut, 'results.mat'));
pl = dlmread(fullfile(dirOut, 'plots.csv'), ',', 1, 0);
dt = dlmread(fullfile(dirOut, 'detections.csv'), ',', 1, 0);
nPl = numel(rn.scans(1).plots.azDeg) + numel(rn.scans(2).plots.azDeg);
nDt = numel(rn.scans(1).det.pulse) + numel(rn.scans(2).det.pulse);
okSave = isequal(sv.result.timeline.startS, rn.timeline.startS) && ...
         isequal(sv.result.scans(2).det.cell, rn.scans(2).det.cell) && ...
         size(pl, 1) == nPl && size(dt, 1) == nDt && ...
         (nDt == 0 || max(abs(dt(:, 5) - [rn.scans(1).det.rangeM; rn.scans(2).det.rangeM])) < 1e-3) && ...
         exist(fullfile(dirOut, 'timeline.csv'), 'file') && exist(fullfile(dirOut, 'summary.txt'), 'file');
ok = check(ok, sprintf('save: results.mat, %d plots and %d detections read back from CSV', ...
           nPl, nDt), okSave);
rmdir(dirOut, 's');

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
