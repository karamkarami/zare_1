function truth = rsp_truth(P, S, nScans)
%RSP_TRUTH Where and how each simulated target should appear.
%
%   truth = rsp_truth(P, S)          first revolution
%   truth = rsp_truth(P, S, nScans)  every revolution 1..nScans
%
%   One entry per target and revolution, when the beam crosses the target:
%   .target .scan                 target number, revolution number
%   .azDeg  .rangeM  .cell        position (range at the crossing time)
%   .velocityMps                  true radial velocity (positive = approaching)
%   .foldedMps                    velocity the radar measures (folded into
%                                 +-lambda*PRF/4)
%   .bin                          Doppler bin of that velocity
%   .pulse                        crossing pulse number
%   .powerDb, .snrDb              echo power and SNR per sample at the beam peak

if nargin < 3 || isempty(nScans)
    nScans = 1;
end
S    = rsp_levels(S);
g    = rsp_geometry(P, S.nRange);
nfft = max(P.fft.nfft, P.fft.nPulses);
perScan = 360 / g.degPerPulse;                                  % pulses per revolution
truth = struct('target', {}, 'scan', {}, 'azDeg', {}, 'rangeM', {}, 'cell', {}, ...
               'velocityMps', {}, 'foldedMps', {}, 'bin', {}, 'pulse', {}, ...
               'powerDb', {}, 'snrDb', {});
for sc = 1:nScans
    for k = 1:numel(S.targets)
        tg = S.targets(k);
        p  = 1 + mod(tg.azDeg - P.radar.azStartDeg, 360) / g.degPerPulse + (sc - 1) * perScan;
        t  = (p - 1) / P.radar.prfHz;
        Rm = tg.rangeM - tg.velocityMps * t;
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
        truth(i).cell        = round((Rm - P.radar.rangeOffsetM) / g.cellM) + 1;
        truth(i).velocityMps = tg.velocityMps;
        truth(i).foldedMps   = (mod(fd + 0.5, 1) - 0.5) * P.radar.prfHz * g.lambdaM / 2;
        truth(i).bin         = b + 1;
        truth(i).pulse       = p;
        truth(i).powerDb     = tg.powerDb;
        truth(i).snrDb       = tg.powerDb - S.noiseDb;
    end
end
end
