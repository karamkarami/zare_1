function truth = rsp_truth(P, S)
%RSP_TRUTH Where and how each simulated target should appear.
%
%   truth = rsp_truth(P, S)
%
%   For every target of S.targets, at the first time the beam crosses it:
%   .azDeg  .rangeM  .cell        position
%   .velocityMps                  true radial velocity (positive = approaching)
%   .foldedMps                    velocity the radar measures (folded into
%                                 +-lambda*PRF/4)
%   .bin                          Doppler bin of that velocity
%   .pulse                        crossing pulse number
%   .snrDb                        SNR per sample at the beam peak

g = rsp_geometry(P, S.nRange);
truth = struct('azDeg', {}, 'rangeM', {}, 'cell', {}, 'velocityMps', {}, ...
               'foldedMps', {}, 'bin', {}, 'pulse', {}, 'snrDb', {});
for k = 1:numel(S.targets)
    tg = S.targets(k);
    p  = 1 + mod(tg.azDeg - P.radar.azStartDeg, 360) / g.degPerPulse;   % crossing pulse
    t  = (p - 1) / P.radar.prfHz;
    Rm = tg.rangeM - tg.velocityMps * t;
    fd = 2 * tg.velocityMps / g.lambdaM / P.radar.prfHz;               % cycles / pulse
    nfft = max(P.fft.nfft, P.fft.nPulses);
    b  = mod(round(fd * nfft), nfft);
    if P.fft.shift
        b = mod(b + floor(nfft/2), nfft);
    end
    truth(k).azDeg       = tg.azDeg;
    truth(k).rangeM      = Rm;
    truth(k).cell        = round((Rm - P.radar.rangeOffsetM) / g.cellM) + 1;
    truth(k).velocityMps = tg.velocityMps;
    truth(k).foldedMps   = (mod(fd + 0.5, 1) - 0.5) * P.radar.prfHz * g.lambdaM / 2;
    truth(k).bin         = b + 1;
    truth(k).pulse       = p;
    truth(k).snrDb       = tg.snrDb;
end
end
