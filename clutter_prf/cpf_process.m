function [y, st, info] = cpf_process(st, blk)
%CPF_PROCESS Clutter-map MTI of one block of pulses; learns the map on the way.
%
%   [y, st, info] = cpf_process(st, blk)
%
%   blk : see cpf_config (x, prfHz, azDeg, rangeM, optional velEN, posEN)
%   y   : filtered pulses x range (same size as blk.x), unit noise gain, ready
%         for the Doppler FFT. Row m uses pulses m, m-1, ... (the previous
%         block's last pulses continue the filter across blocks).
%   info: .filter  bank index of every sample (0 = not filtered)
%         .cell    map cell of every sample (0 = outside the map)
%         .noise   noise power in use
%         .scanDone true when the antenna passed north in this block
%
%   Steps
%   1. Platform Doppler of the ground, 2 v.u / lambda, is removed pulse by
%      pulse (cfg.mti.motionComp), so fixed ground clutter is at 0 Hz.
%   2. Every sample is filtered with the bank filter of its map cell (clutter
%      power, width and mean velocity learned on earlier scans, plus the scan
%      and platform spread now), shifted to the cell's mean velocity.
%      Cells with fewer than cfg.map.minUpdates updates use cfg.clutter.
%   3. Correlations at lags 0 ... cfg.map.maxLag are accumulated per map cell.
%      A cell is folded into the map once the beam has left it.

cfg = st.cfg;
lam = st.lambda;
x = blk.x;
[M, Nc] = size(x);
prf = blk.prfHz;
T = 1 / prf;
az = mod(blk.azDeg(:), 360);
rng = blk.rangeM(:).';
vel = getRows(blk, 'velEN', M, [0 0]);
if isfield(blk, 'posEN') && ~isempty(blk.posEN)
    pos = getRows(blk, 'posEN', M, [0 0]);
else
    pos = st.pos + [0 0; cumsum(vel(1:end-1, :) * T, 1)];
end
st.pos = pos(end, :) + vel(end, :) * T;
st.vel = vel(end, :);

% pending cells were measured at another PRF: fold them in first
if ~isempty(st.accPrf) && abs(st.accPrf - prf) > 1e-9
    st = cpf_map_commit(st, find(st.acc.n0 > 0), 1 / st.accPrf);
    st.tail = [];                                % lags must not span two PRFs
end
st.accPrf = prf;
if isempty(st.map.origin)
    st.map.origin = pos(1, :);
elseif norm(pos(end, :) - st.map.origin) > cfg.map.recenterM
    st = cpf_map_commit(st, find(st.acc.n0 > 0), T);
    st = cpf_map_recenter(st, pos(end, :));
    st.tail = [];
end

% antenna speed
if M > 1
    d = abs(diff(unwrap(az * pi / 180))) * 180 / pi;
    om = median(d) * prf;
    if om > 0
        st.omegaDps = om;
    end
end
[sScan, kMot] = cpf_spread(st);

% 1. platform Doppler
u = [sind(az) cosd(az)];
vr    = sum(vel .* u, 2);                          % platform speed towards the beam
vPerp = vel(:, 1) .* cosd(az) - vel(:, 2) .* sind(az);
sMot2 = (kMot * vPerp).^2;
if cfg.mti.motionComp
    phi = st.phi + cumsum(2 * pi * 2 * vr / lam * T);
else
    phi = zeros(M, 1);
end
st.phi = mod(phi(end), 2 * pi);
xd = x .* exp(-1j * phi);

% noise
if isempty(cfg.noise.power)
    pc = sort(mean(abs(x).^2, 1));
    low = pc(max(1, round(cfg.noise.fraction * Nc)));
    nEst = mean(pc(pc <= 2 * low));             % all cells near the weakest ones
    if isempty(st.noise)
        st.noise = nEst;
    else
        a = 1 - (1 - cfg.noise.alpha)^(M / 256);  % weight by the block length
        st.noise = (1 - a) * st.noise + a * nEst;
    end
else
    st.noise = cfg.noise.power;
end

L = cpf_map_index(st, pos, az, rng);

% 2. filter
nT = cfg.mti.maxLength;
K = cfg.map.maxLag;
nTail = max(nT - 1, K);
if isempty(st.tail) || size(st.tail.x, 2) ~= Nc
    st.tail.x = zeros(0, Nc);
    st.tail.L = zeros(0, Nc);
end
xe = [st.tail.x; xd];
Le = [st.tail.L; L];
p0 = size(st.tail.x, 1);
[b, q] = cellFilter(st, L, sMot2, sScan, T);
[bank, st] = getBank(st, prf, b);
y = complex(zeros(M, Nc));
nb = size(bank.W, 2) + 1;
key = b + nb * (q - min(q(:)));
[keys, ~, ik] = unique(key(:));
ik = reshape(ik, M, Nc);
k = (0:nT-1).';
for j = 1:numel(keys)
    mask = ik == j;
    cols = any(mask, 1);
    bj = b(find(mask, 1));
    if bj == 0
        y(mask) = xd(mask);
        continue
    end
    qj = q(find(mask, 1));
    w = bank.W(:, bj) .* exp(2j * pi * (2 * qj * cfg.mti.vmStepMps / lam) * T * k);
    ye = filter(w, 1, xe(:, cols), [], 1);
    yj = complex(zeros(M, Nc));
    yj(:, cols) = ye(p0+1:end, :);
    y(mask) = yj(mask);
end
if p0 < nT - 1
    y(1:min(nT - 1 - p0, M), :) = 0;             % start of the data: filter not filled yet
end
if strcmpi(cfg.mti.output, 'radar')
    y = y .* exp(1j * phi);
end

% 3. accumulate correlations
v0 = L > 0;
i0 = L(v0);
st.acc.s0 = st.acc.s0 + accumarray(i0(:), abs(xd(v0(:))).^2, size(st.acc.s0));
st.acc.n0 = st.acc.n0 + accumarray(i0(:), 1, size(st.acc.n0));
mot = repmat(sMot2, 1, Nc);
st.acc.m2 = st.acc.m2 + accumarray(i0(:), mot(v0(:)), size(st.acc.m2));
for k = 1:K
    r = (1:M) + p0 - k;                          % rows of xe k pulses earlier
    ok = r >= 1;
    if ~any(ok)
        break
    end
    vk = false(M, Nc);
    vk(ok, :) = v0(ok, :) & Le(r(ok), :) > 0;
    prev = complex(zeros(M, Nc));
    prev(ok, :) = xe(r(ok), :);
    ik = L(vk);
    rk = xd(vk(:)) .* conj(prev(vk(:)));
    st.acc.s(:, k) = st.acc.s(:, k) + accumarray(ik(:), real(rk), [size(st.acc.s, 1) 1]) ...
                     + 1j * accumarray(ik(:), imag(rk), [size(st.acc.s, 1) 1]);
    st.acc.n(:, k) = st.acc.n(:, k) + accumarray(ik(:), 1, [size(st.acc.s, 1) 1]);
end

% fold in the cells the beam has left (not touched by the last azimuth cell of pulses)
dAz = abs(mod(az - az(end) + 180, 360) - 180);
held = L(dAz < cfg.map.azCellDeg, :);
busy = false(size(st.acc.n0));
busy(held(held > 0)) = true;
st = cpf_map_commit(st, find(st.acc.n0 > 0 & ~busy), T);

% tail for the next block
keep = min(nTail, size(xe, 1));
st.tail.x = xe(end-keep+1:end, :);
st.tail.L = Le(end-keep+1:end, :);

% revolutions
azAll = [st.lastAz; az];
wraps = sum(diff(azAll) < -180);
st.scan = st.scan + wraps;
st.lastAz = az(end);

info.filter = b;
info.cell = L;
info.noise = st.noise;
info.scanDone = wraps > 0;
end

function v = getRows(blk, name, M, def)
if isfield(blk, name) && ~isempty(blk.(name))
    v = blk.(name);
else
    v = def;
end
if size(v, 1) == 1
    v = repmat(v, M, 1);
end
end

function [bank, st] = getBank(st, prf, b)
% filter bank of this PRF (cached), with the filters needed by b designed
i = find(abs(st.bankPrf - prf) < 1e-6, 1);
if isempty(i)
    st.banks{end+1} = cpf_bank(st.cfg, prf, st.lambda);
    st.bankPrf(end+1) = prf;
    keepN = 3 * st.cfg.prf.count;               % cache size
    if numel(st.banks) > keepN
        st.banks = st.banks(end-keepN+1:end);
        st.bankPrf = st.bankPrf(end-keepN+1:end);
    end
    i = numel(st.banks);
end
st.banks{i} = cpf_bank(st.cfg, st.banks{i}, st.lambda, b(:).');
bank = st.banks{i};
end

function [b, q] = cellFilter(st, L, sMot2, sScan, T)
% bank index (0 = no filter) and mean velocity step of every sample
cfg = st.cfg;
m = st.map;
cnrDb = cfg.clutter.cnrDb * ones(size(L));
sig2  = cfg.clutter.sigmaVMps^2 * ones(size(L));
vm    = zeros(size(L));
in = L > 0;
own = false(size(L));
own(in) = m.n(L(in)) >= cfg.map.minUpdates;
if ~isempty(st.noise) && st.noise > 0
    cnrDb(own) = 10 * log10(max(m.C(L(own)), 1e-30) / st.noise);
end
sig2(own) = m.sig2(L(own));
vm(own)   = m.vm(L(own));
sTot = cfg.mti.widthMargin * T * sqrt((2 / st.lambda)^2 * sig2 + sScan^2 + repmat(sMot2, 1, size(L, 2)));
sG = cfg.mti.sGrid;
cG = cfg.mti.cnrGridDb;
iS = gridUp(sG, sTot);
iC = gridUp(cG, cnrDb + cfg.mti.cnrMarginDb);
b = iS + (iC - 1) * numel(sG);
b(cnrDb < cfg.mti.bypassCnrDb) = 0;
q = round(vm / cfg.mti.vmStepMps);
end

function i = gridUp(g, x)
% index of the first grid value >= x (clamped to the grid)
x = min(max(x, g(1)), g(end));
i = ceil(interp1(g, 1:numel(g), x) - 1e-9);
i = min(max(i, 1), numel(g));
end
