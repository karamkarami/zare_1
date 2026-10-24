function [blk, sim, truth] = csim_block(sim, prfHz)
%CSIM_BLOCK Next block of pulses at the given PRF.
%
%   [blk, sim, truth] = csim_block(sim, prfHz)
%
%   blk   : x (pulses x range, complex), prfHz, azDeg, rangeM, velEN, timeS
%   truth : .scanDone  the block ends the revolution (the antenna reached north)
%           .scan      revolution of this block (1, 2, ...)
%           .tgt       per target at the block middle: rangeM, azDeg (from the
%                      radar), vRadMps (over ground, approaching > 0), cell
%           .tgtPos    pulses x 2 x nTargets target positions
%   Change sim.vel between blocks to make the platform move.
%   Range response: an echo at range R is shared between the two nearest
%   range cells (triangular response of one cell width), so echoes move
%   smoothly from cell to cell when the radar or the target moves.

S = sim.S;
lam = sim.lambda;
T = 1 / prfHz;
endAz = 360 * (floor((sim.azAbs - S.azStartDeg) / 360 + 1e-9) + 1) + S.azStartDeg;
nLeft = max(ceil((endAz - sim.azAbs) / (sim.omega * T) - 1e-9), 1);
M = min(S.blockPulses, nLeft);
k = (0:M-1).';
t = sim.t + k * T;
azAbs = sim.azAbs + sim.omega * k * T;
az = mod(azAbs, 360);
pos = sim.pos + (t - sim.t) * sim.vel;
rangeM = (1:S.nRange) * S.rangeCellM;
x = complex(zeros(M, S.nRange));
bw = S.beamwidthDeg;

% clutter, 32 pulses at a time
for i0 = 1:32:M
    ii = i0:min(i0 + 31, M);
    pm = pos(ii(ceil(end / 2)), :);
    dE = sim.cl.E - pm(1);  dN = sim.cl.N - pm(2);
    aRel = atan2(dE, dN) * 180 / pi;
    span = sim.omega * T * numel(ii);
    sel = find(abs(mod(aRel - az(ii(ceil(end / 2))) + 180, 360) - 180) < 2 * bw + span);
    if isempty(sel)
        continue
    end
    E = sim.cl.E(sel).' - pos(ii, 1);              % pulses x scatterers
    N = sim.cl.N(sel).' - pos(ii, 2);
    R = hypot(E, N);
    th = mod(atan2(E, N) * 180 / pi - az(ii) + 180, 360) - 180;
    g = exp(-4 * log(2) * th.^2 / bw^2);          % two-way voltage
    nC = size(sim.cl.f, 2);
    fl = zeros(numel(ii), numel(sel));
    for c = 1:nC
        fl = fl + sim.cl.c(sel, c).' .* exp(2j * pi * t(ii) * sim.cl.f(sel, c).');
    end
    a = (sim.cl.amp(sel).' .* g) .* fl .* exp(-4j * pi * R / lam);
    r0 = floor(R(ceil(end / 2), :) / S.rangeCellM);
    fr = min(max(R / S.rangeCellM - r0, 0), 1);     % share of the upper cell
    for d = 0:1
        rc = r0 + d;
        ok = rc >= 1 & rc <= S.nRange;
        Sp = sparse(find(ok), rc(ok), 1, numel(sel), S.nRange);
        if d == 0
            x(ii, :) = x(ii, :) + (a .* (1 - fr)) * Sp;
        else
            x(ii, :) = x(ii, :) + (a .* fr) * Sp;
        end
    end
end

% targets
nT = size(sim.tg.pos, 1);
truth.tgtPos = zeros(M, 2, nT);
for j = 1:nT
    p = sim.tg.pos(j, :) + (t - sim.t) * sim.tg.vel(j, :);
    truth.tgtPos(:, :, j) = p;
    E = p(:, 1) - pos(:, 1);  N = p(:, 2) - pos(:, 2);
    R = hypot(E, N);
    th = mod(atan2(E, N) * 180 / pi - az + 180, 360) - 180;
    g = exp(-4 * log(2) * th.^2 / bw^2);
    a = sqrt(S.noisePower * 10^(sim.tg.snrDb(j) / 10)) * sim.tg.amp(j) * g .* exp(-4j * pi * R / lam);
    r0 = floor(R / S.rangeCellM);
    fr = R / S.rangeCellM - r0;
    for d = 0:1
        rc = r0 + d;
        ok = rc >= 1 & rc <= S.nRange & g > 1e-4;
        idx = find(ok);
        sh = (1 - fr(idx)) * (d == 0) + fr(idx) * (d == 1);
        x(sub2ind(size(x), idx, rc(idx))) = x(sub2ind(size(x), idx, rc(idx))) + a(idx) .* sh;
    end
    mid = ceil(M / 2);
    u = [E(mid) N(mid)] / R(mid);
    truth.tgt(j).rangeM  = R(mid);
    truth.tgt(j).azDeg   = mod(atan2(E(mid), N(mid)) * 180 / pi, 360);
    truth.tgt(j).vRadMps = -sim.tg.vel(j, :) * u.';
    truth.tgt(j).cell    = round(R(mid) / S.rangeCellM);
end

x = x + sqrt(S.noisePower / 2) * (randn(M, S.nRange) + 1j * randn(M, S.nRange));

blk.x = x;
blk.prfHz = prfHz;
blk.azDeg = az;
blk.rangeM = rangeM;
blk.velEN = sim.vel;
blk.timeS = t;

truth.scan = sim.scan + 1;
truth.scanDone = M == nLeft;
sim.t = t(end) + T;
sim.azAbs = azAbs(end) + sim.omega * T;
sim.tg.pos = sim.tg.pos + (sim.t - t(1)) * sim.tg.vel;
sim.pos = pos(end, :) + sim.vel * T;
if truth.scanDone
    sim.scan = sim.scan + 1;
    if S.swerling
        sim.tg.amp = (randn(nT, 1) + 1j * randn(nT, 1)) / sqrt(2);
    end
end
end
