function truth = rsp_truth(P, S, nScans)
%RSP_TRUTH Where and how each simulated target should appear.
%
%   truth = rsp_truth(P, S)          first revolution
%   truth = rsp_truth(P, S, nScans)  every revolution 1..nScans
%
%   One entry per target and revolution, when the beam crosses the target:
%   .target .scan                 target number, revolution number
%   .azDeg  .rangeM  .cell        position (range at the crossing time)
%   .apparentRangeM               range where the echo appears: beyond the
%                                 unambiguous range c/2PRF it folds back
%                                 (.cell is the apparent cell)
%   .ambiguous                    true = second-time-around echo
%   .timeS                        crossing time (P.radar.startTimeS included)
%   .velocityMps                  true radial velocity (positive = approaching)
%   .foldedMps                    velocity the radar measures (folded into
%                                 +-lambda*PRF/4)
%   .bin                          Doppler bin of that velocity
%   .pulse                        crossing pulse number
%   .powerDb, .snrDb              echo power and SNR per sample at the beam peak
%                                 (the SNR includes the terrain shadow)
%   .shadowDb                     loss behind a mountain layer (0 = visible)

if nargin < 3 || isempty(nScans)
    nScans = 1;
end
S    = rsp_levels(S);
g    = rsp_geometry(P, S.nRange);
nfft = max(P.fft.nfft, P.fft.nPulses);
perScan = 360 / g.degPerPulse;                                  % pulses per revolution
truth = struct('target', {}, 'scan', {}, 'azDeg', {}, 'rangeM', {}, 'cell', {}, ...
               'apparentRangeM', {}, 'ambiguous', {}, 'timeS', {}, ...
               'velocityMps', {}, 'foldedMps', {}, 'bin', {}, 'pulse', {}, ...
               'powerDb', {}, 'snrDb', {}, 'shadowDb', {});
L = [];
if isfield(S, 'clutter') && ~isempty(S.clutter)
    L = rsp_clutter_types(S.clutter);
end
t0 = 0;
if isfield(P.radar, 'startTimeS') && ~isempty(P.radar.startTimeS)
    t0 = P.radar.startTimeS;
end
for sc = 1:nScans
    for k = 1:numel(S.targets)
        tg = S.targets(k);
        p  = 1 + mod(tg.azDeg - P.radar.azStartDeg, 360) / g.degPerPulse + (sc - 1) * perScan;
        t  = t0 + (p - 1) / P.radar.prfHz;
        Rm = tg.rangeM - tg.velocityMps * t;
        Ra = P.radar.rangeOffsetM + mod(Rm - P.radar.rangeOffsetM, g.rangeUnambM);
        sh = 0;
        if ~isempty(L)
            sh = 10*log10(rsp_clutter_shadow(L, tg.azDeg, Rm));
        end
        fd = 2 * tg.velocityMps / g.lambdaM / P.radar.prfHz;           % cycles / pulse
        b  = mod(round(fd * nfft), nfft);
        if P.fft.shift
            b = mod(b + floor(nfft/2), nfft);
        end
        i = numel(truth) + 1;
        truth(i).target      = k;
        truth(i).scan        = sc;
        truth(i).azDeg       = tg.azDeg;
        truth(i).rangeM      = Rm;
        truth(i).cell        = round((Ra - P.radar.rangeOffsetM) / g.cellM) + 1;
        truth(i).apparentRangeM = Ra;
        truth(i).ambiguous   = abs(Ra - Rm) > 1;
        truth(i).timeS       = t;
        truth(i).velocityMps = tg.velocityMps;
        truth(i).foldedMps   = (mod(fd + 0.5, 1) - 0.5) * P.radar.prfHz * g.lambdaM / 2;
        truth(i).bin         = b + 1;
        truth(i).pulse       = p;
        truth(i).powerDb     = tg.powerDb;
        truth(i).snrDb       = tg.powerDb - S.noiseDb + sh;
        truth(i).shadowDb    = sh;
    end
end
end
