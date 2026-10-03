function G = rsp_antenna_pattern(dAzDeg, A)
%RSP_ANTENNA_PATTERN One-way azimuth power pattern of the antenna (linear, peak 1).
%
%   G = rsp_antenna_pattern(dAzDeg, A)
%
%   dAzDeg : angle from the beam axis [deg], any size, any range (wrapped)
%   A      : P.antenna
%            .beamwidthDeg  3 dB beamwidth                       (4)
%            .sidelobeDb    first sidelobe level, one way          (-25)
%            .backlobeDb    level from 90 deg to the back, one way (-40)
%            .nullDepthDb   depth of the nulls between sidelobes   (-20)
%
%   Main lobe : Gaussian with the given 3 dB width.
%   Sidelobes : lobes one beamwidth apart; the envelope starts at sidelobeDb
%               on the first sidelobe (1.6 beamwidths) and falls (linear in
%               log-angle) to backlobeDb at 90 deg, constant behind.
%   The two-way echo power of a point target is G^2 (amplitude G).

th  = abs(mod(dAzDeg + 180, 360) - 180);              % 0 .. 180 deg
bw  = A.beamwidthDeg;
gm  = exp(-4*log(2) * (th/bw).^2);                    % main lobe

th1 = 1.6 * bw;                                       % first sidelobe
x   = min(max(log(max(th, th1)/th1) / log(90/th1), 0), 1);
env = A.sidelobeDb + (A.backlobeDb - A.sidelobeDb) * x;          % dB
rip = cos(pi*(th - th1)/bw).^2;                                  % lobes and nulls
rip = max(rip, 10^(A.nullDepthDb/10));
gs  = 10.^(env/10) .* rip;
gs(th < th1 - bw/2) = 0;                              % inside the main lobe

G = max(gm, gs);
end
