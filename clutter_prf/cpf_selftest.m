function ok = cpf_selftest()
%CPF_SELFTEST Checks of the clutter-map MTI and the PRF selection.
%
%   ok = cpf_selftest()
%
%   1. filter design : unit norm, clutter residue under the limit, notch at the
%                      clutter, a shifted filter has its notch at the shift
%   2. estimation    : power, intrinsic spread and mean velocity of synthetic
%                      Gaussian clutter come back from the map
%   3. motion        : a moving platform's clutter is cancelled with motion
%                      compensation and not without it
%   4. blocks        : the output does not depend on how the pulses are split
%                      into blocks
%   5. PRF selection : local search = exhaustive search; limits and cost
%   6. map           : re-centring keeps what was learned; PRFs used in turn
%
%   Uses synthetic data only (independent of clutter_sim).

rng(3);
ok = true;
fprintf('cpf_selftest\n');
cfg = cpf_config();
lam = 299792458 / cfg.radar.fcHz;

%% 1. Filter design ---------------------------------------------------------------------
prf = 1000;
B = cpf_band(cfg, prf, lam);
s = 3 / prf;
[w, caDb] = cpf_design_filter(B, s, 60, cfg);
rc = exp(-2 * pi^2 * s^2 * (0:numel(w)-1).^2);
res = 10 * log10(10^6 * real(w' * toeplitz(rc) * w));
ok = check(ok, sprintf('design: unit norm, residue %.1f dB <= %g dB, CA %.1f dB', res, ...
           cfg.mti.residueDb, caDb), abs(norm(w) - 1) < 1e-12 && res <= cfg.mti.residueDb + 1e-6);
g0 = abs(sum(w))^2;
gp = mean(abs(B.S * w).^2);
ok = check(ok, sprintf('design: gain at 0 Hz %.0f dB, mean target gain %.1f dB', ...
           10 * log10(g0), 10 * log10(gp)), 10 * log10(g0) < -60 && 10 * log10(gp) > -3);
f0 = 37;
ws = w .* exp(2j * pi * f0 / prf * (0:numel(w)-1).');
ok = check(ok, 'design: shifted taps put the notch at the shift', ...
           abs(cpf_steering(f0, prf, numel(w)) * ws)^2 < 1e-6);
cfgA = cfg;  cfgA.mti.lengthMode = 'auto';
wA = cpf_design_filter(B, s, 20, cfgA);
ok = check(ok, sprintf('design: auto length uses %d taps for CNR 20 dB', nnz(abs(wA) > 0)), ...
           nnz(abs(wA) > 0) < cfg.mti.maxLength);

%% 2. Estimation ------------------------------------------------------------------------
c2 = cfg;
c2.radar.rpm = 0;  c2.noise.power = 1;  c2.map.minUpdates = 1;
st = cpf_init(c2, [900 1000 1100 1200]);
M = 4000;  nC = 64;
sigV = 0.6;  vm = 1.5;  cnrDb = 40;
X = gaussClutter(M, nC, 2 * sigV / lam, 2 * vm / lam, prf, cnrDb) + noise(M, nC);
blk = struct('x', X, 'prfHz', prf, 'azDeg', 45 * ones(M, 1), 'rangeM', 2000 + (0:nC-1) * 30);
[~, st] = cpf_process(st, blk);
st = cpf_map_commit(st, find(st.acc.n0 > 0), 1 / prf);
cells = find(st.map.n > 0);
Cdb = 10 * log10(mean(st.map.C(cells)));
sv = sqrt(mean(st.map.sig2(cells)));
v  = mean(st.map.vm(cells));
ok = check(ok, sprintf('estimate: CNR %.2f dB (40), sigma_v %.3f m/s (%.1f), mean %.3f m/s (%.1f)', ...
           Cdb, sv, sigV, v, vm), abs(Cdb - cnrDb) < 0.5 && abs(sv / sigV - 1) < 0.1 && abs(v - vm) < 0.05);

%% 3. Platform motion -------------------------------------------------------------------
% map learned while stationary, then the platform moves: the ground clutter
% moves to the platform Doppler
velEN = [12 -5];  azDeg = 60;
fp = 2 * (velEN * [sind(azDeg); cosd(azDeg)]) / lam;
h = M / 2;
X1 = gaussClutter(h, nC, 2 * 0.2 / lam, 0, prf, 50) + noise(h, nC);
X2 = gaussClutter(h, nC, 2 * 0.2 / lam, fp, prf, 50) + noise(h, nC);
b1 = struct('x', X1, 'prfHz', prf, 'azDeg', azDeg * ones(h, 1), 'rangeM', 2000 + (0:nC-1) * 30);
b2 = b1;  b2.x = X2;  b2.velEN = velEN;
r = zeros(1, 2);
for comp = [true false]
    c3 = c2;  c3.mti.motionComp = comp;
    st = cpf_init(c3, [900 1000 1100 1200]);
    [~, st] = cpf_process(st, b1);
    st = cpf_map_commit(st, find(st.acc.n0 > 0), 1 / prf);
    [y, st] = cpf_process(st, b2);
    r(2 - comp) = 10 * log10(mean(mean(abs(y(cfg.mti.maxLength:end, :)).^2)));
end
ok = check(ok, sprintf('motion: %.0f Hz ground Doppler, residue %.1f dB with compensation, %.1f dB without', ...
           fp, r(1), r(2)), r(1) < 1 && r(2) > 20);

%% 4. Block size ------------------------------------------------------------------------
c4 = cfg;  c4.noise.power = 1;
st0 = cpf_init(c4, [900 1000 1100 1200]);
M = 3000;  nC = 40;
az = mod(10 + (0:2*M-1).' * 0.036, 360);
X = gaussClutter(2 * M, nC, 2 * 0.3 / lam, 0, prf, 45) + noise(2 * M, nC);
blk = struct('x', X(1:M, :), 'prfHz', prf, 'azDeg', az(1:M), 'rangeM', 1000 + (0:nC-1) * 150);
st = st0;  [~, st] = cpf_process(st, blk);              % learn a map first
X = X(M+1:end, :);  az = az(M+1:end);                   % then go on turning
blk.x = X;  blk.azDeg = az;
stA = st;  stB = st;
[yA, stA] = cpf_process(stA, blk);
cuts = [0 1 7 300 301 1700 M];
yB = zeros(M, nC);
for i = 1:numel(cuts) - 1
    rr = cuts(i)+1:cuts(i+1);
    b = blk;  b.x = X(rr, :);  b.azDeg = az(rr);
    [yB(rr, :), stB] = cpf_process(stB, b);
end
ok = check(ok, 'blocks: one block = 6 blocks of 1 ... 1399 pulses (output and map)', ...
           max(abs(yA(:) - yB(:))) < 1e-9 * max(abs(yA(:))) && ...
           max(abs(stA.map.C - stB.map.C)) < 1e-9 * max(stA.map.C));

%% 5. PRF selection ---------------------------------------------------------------------
c5 = cfg;  c5.prf.stepHz = 50;  c5.prf.count = 3;
d.cnrDb = 60;  d.sigmaHz = 4;  d.speedMps = 0;
e = cpf_select_prfs(c5, d, lam);
c5.prf.maxExhaustive = 0;
l = cpf_select_prfs(c5, d, lam);
ok = check(ok, sprintf('PRF: local search %s = exhaustive %s', mat2str(l.prfHz), mat2str(e.prfHz)), ...
           abs(l.J - e.J) < 1e-9);
ev = cpf_prf_eval(c5, d, e.prfHz, lam);
ok = check(ok, 'PRF: cpf_prf_eval of the chosen set gives the same cost', abs(ev.J - e.J) < 1e-9);
bl = sum(e.blind, 1);
ok = check(ok, sprintf('PRF: blind PRFs at once %d (limit %d), CA %.1f >= %.1f dB', e.blindMax, ...
           c5.prf.maxBlind, e.caWorstDb, e.caReqDb), e.feasible && max(bl) == e.blindMax);
same = cpf_prf_eval(c5, d, [1000 1000 1000], lam);
ok = check(ok, sprintf('PRF: three equal PRFs are blind together (%d)', same.blindMax), same.blindMax == 3);

%% 6. Map and PRF order -----------------------------------------------------------------
st = st0;  st.map.origin = [0 0];
st.map.C(:) = 1:numel(st.map.C);
p = [3000 4000];
L0 = cpf_map_index(st, [0 0], atan2d(p(1), p(2)), norm(p));
st2 = cpf_map_recenter(st, [400 -300]);
L1 = cpf_map_index(st2, [400 -300], atan2d(p(1) - 400, p(2) + 300), norm(p - [400 -300]));
ok = check(ok, 're-centre: a ground point keeps its clutter', st2.map.C(L1) == st.map.C(L0));
st = st0;  q = zeros(1, 2 * numel(st.prfSet));
for i = 1:numel(q)
    [q(i), st] = cpf_next_prf(st, 0);
end
ok = check(ok, 'cpf_next_prf: PRFs used in turn', isequal(q, [st.prfSet st.prfSet]));

if ok
    fprintf('all checks passed\n');
else
    fprintf('SOME CHECKS FAILED\n');
end
end

function X = gaussClutter(M, nC, sigmaHz, f0, prf, cnrDb)
% Gaussian spectrum clutter, independent per range cell, power 10^(cnrDb/10)
f = ((0:M-1).' - floor(M / 2)) * prf / M;
f = ifftshift(f);
F = mod(f - f0 + prf / 2, prf) - prf / 2;
H = exp(-F.^2 / (4 * sigmaHz^2));
H = H / sqrt(mean(H.^2));
X = ifft(fft(noise(M, nC)) .* H) * sqrt(10^(cnrDb / 10));
end

function n = noise(M, nC)
n = (randn(M, nC) + 1j * randn(M, nC)) / sqrt(2);
end

function ok = check(ok, name, pass)
if pass
    fprintf('  ok    %s\n', name);
else
    fprintf('  FAIL  %s\n', name);
    ok = false;
end
end
