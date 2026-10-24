function st = cpf_map_commit(st, cells, T)
%CPF_MAP_COMMIT Fold the accumulated correlations of some cells into the map.
%
%   st = cpf_map_commit(st, cells, T)      T = pulse interval of the data
%
%   Per cell: clutter power C = R0 - noise. Above cfg.map.minCnrDbWidth also
%   the spectrum. For a Gaussian spectrum |R(k)| = C exp(-2 pi^2 sigma^2 T^2 k^2),
%   so sigma comes from a straight-line fit of log|R(k)| against k^2 over the
%   lags 1 ... cfg.map.maxLag that are still well correlated (narrow clutter
%   hardly decorrelates in one pulse, so one lag alone is far too noisy).
%   Total sigma minus the scan and platform spread gives the intrinsic sigma_v;
%   the phase of R(1) gives the mean velocity. Both are in m/s, so a cell
%   learned at one PRF serves every PRF.
%   First update: taken as measured. Later: recursive mean with weight
%   cfg.map.alpha, power in dB and limited to +-cfg.map.clipDb per update,
%   spread limited to a factor 4 per update (a target passing through a cell
%   does not wipe out what was learned).

if isempty(cells)
    return
end
cfg = st.cfg;
lam = st.lambda;
a = st.acc;
m = st.map;
c = cells(:);
N = st.noise;
R0 = a.s0(c) ./ a.n0(c);
Rk = a.s(c, :) ./ max(a.n(c, :), 1);
R1 = Rk(:, 1);
C  = max(R0 - N, 1e-6 * N);
[sScan, ~] = cpf_spread(st);
mot2 = a.m2(c) ./ a.n0(c);
ok = C >= N * 10^(cfg.map.minCnrDbWidth / 10) & a.n(c, 1) > 0;
sf2 = widthFit(Rk, a.n(c, :), C, T);                    % total sigma^2 [Hz^2]
ok = ok & isfinite(sf2);
sig2 = max((lam / 2)^2 * (sf2 - sScan^2 - mot2), cfg.clutter.sigmaVMinMps^2);
vm = lam / 2 * angle(R1) / (2 * pi * T);

al = cfg.map.alpha;
first = m.n(c) == 0;
dB = 10 * log10(C ./ max(m.C(c), 1e-30));
dB = min(max(dB, -cfg.map.clipDb), cfg.map.clipDb);
Cn = m.C(c) .* 10.^(al * dB / 10);
Cn(first) = C(first);
m.C(c) = Cn;

w1 = ok & m.nw(c) == 0;
wn = ok & m.nw(c) > 0;
s2 = m.sig2(c);
v  = m.vm(c);
s2(w1) = sig2(w1);
v(w1)  = vm(w1);
lim = min(max(sig2(wn), s2(wn) / 4), 4 * s2(wn));
s2(wn) = (1 - al) * s2(wn) + al * lim;
v(wn)  = (1 - al) * v(wn) + al * vm(wn);
m.sig2(c) = s2;
m.vm(c) = v;
m.n(c)  = m.n(c) + 1;
m.nw(c) = m.nw(c) + ok;
st.map = m;

st.acc.s0(c) = 0;  st.acc.n0(c) = 0;
st.acc.s(c, :) = 0;  st.acc.n(c, :) = 0;
st.acc.m2(c) = 0;
end

function sf2 = widthFit(Rk, nk, C, T)
% Weighted line fit of log|R(k)| = a - 2 pi^2 sigma^2 T^2 k^2 per cell (row).
% Lags with |R(k)| below 0.1 C (decorrelated) or without samples are left out;
% with fewer than 3 lags left, lag 1 alone.
K = size(Rk, 2);
x = (1:K).^2;
y = log(max(abs(Rk), 1e-300));
w = nk .* (abs(Rk) > 0.1 * C);
sw = sum(w, 2);
xm = (w * x.') ./ max(sw, 1);
ym = sum(w .* y, 2) ./ max(sw, 1);
dx = x - xm;
b = -sum(w .* dx .* (y - ym), 2) ./ max(sum(w .* dx.^2, 2), 1e-30);
sf2 = max(b, 0) / (2 * pi^2 * T^2);
few = sum(w > 0, 2) < 3;
rho = min(abs(Rk(:, 1)) ./ C, 1 - 1e-12);
sf2(few) = -log(max(rho(few), 1e-12)) / (2 * pi^2 * T^2);
sf2(sum(nk > 0, 2) == 0) = NaN;
end
