function sim = csim_init(S)
%CSIM_INIT Draw the clutter scatterers and set up targets and antenna.
%
%   sim = csim_init(S)        S = csim_config()

rng(S.seed);
sim.S = S;
sim.lambda = 299792458 / S.fcHz;
sim.t = 0;
sim.azAbs = S.azStartDeg;                 % unwrapped antenna azimuth [deg]
sim.omega = 6 * S.rpm;                    % [deg/s]
sim.pos = [0 0];
sim.vel = S.platform.velEN;
sim.scan = 0;

% clutter scatterers on a polar grid around the start point
C = S.clutter;
r = (ceil(C.minRangeM / S.rangeCellM):floor(C.maxRangeM / S.rangeCellM)) * S.rangeCellM;
a = 0:C.azStepDeg:360 - C.azStepDeg;
[A, R] = ndgrid(a, r);
A = A(:);
R = R(:) + (rand(numel(R), 1) - 0.5) * S.rangeCellM;    % anywhere inside its cell
pDb = C.cnrDb - 10 * C.rangeExp * log10(R / C.refRangeM) + C.textureDb * randn(size(R));
for k = 1:size(C.patches, 1)
    p = C.patches(k, :);
    in = abs(mod(A - p(1) + 180, 360) - 180) <= p(2) / 2 & abs(R - p(3)) <= p(4) / 2;
    pDb(in) = pDb(in) + p(5);
end
beamSum = S.beamwidthDeg * sqrt(pi / (8 * log(2))) / C.azStepDeg;   % sum of g^2 over the scatterers in the beam
sim.cl.E = R .* sind(A);
sim.cl.N = R .* cosd(A);
sim.cl.amp = sqrt(S.noisePower * 10.^(pDb / 10) / beamSum);
nC = C.components;
sim.cl.f = (2 * C.sigmaVMps / sim.lambda) * randn(numel(R), nC);       % [Hz]
sim.cl.c = (randn(numel(R), nC) + 1j * randn(numel(R), nC)) / sqrt(2 * nC);
sim.cl.meanCnrDb = pDb;

% targets
T = S.targets;
sim.tg.pos = zeros(numel(T), 2);
sim.tg.vel = zeros(numel(T), 2);
for k = 1:numel(T)
    u = [sind(T(k).azDeg) cosd(T(k).azDeg)];
    sim.tg.pos(k, :) = T(k).rangeM * u;
    sim.tg.vel(k, :) = -T(k).vRadMps * u + T(k).vTanMps * [u(2) -u(1)];
end
sim.tg.snrDb = [T.snrDb].';
sim.tg.amp = swerling(numel(T), S.swerling);
end

function a = swerling(n, on)
if on
    a = (randn(n, 1) + 1j * randn(n, 1)) / sqrt(2);
else
    a = ones(n, 1);
end
end
