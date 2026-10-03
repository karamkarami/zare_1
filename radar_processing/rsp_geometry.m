function g = rsp_geometry(P, nRange, pulses)
%RSP_GEOMETRY Physical axes of the data: range, azimuth, time, velocity.
%
%   g = rsp_geometry(P, nRange, pulses)
%
%   P       : parameters (P.fs, P.radar, P.fft)
%   nRange  : range cells per pulse
%   pulses  : global pulse numbers (1 = first pulse of the scan), optional
%
%   g.cellM        range cell size [m]
%   g.rangeM       range of every cell [m] (cell 1 = P.radar.rangeOffsetM)
%   g.lambdaM      wavelength [m]
%   g.vUnambMps    unambiguous velocity span, lambda*PRF/2 [m/s]
%   g.binVelMps    radial velocity at the centre of each Doppler bin [m/s]
%                  (positive = approaching), folded into +-vUnamb/2
%   g.degPerPulse  antenna rotation between pulses [deg]
%   g.outputDelayPulses  lag of the centre of the data behind an output row
%                  of the chain (plots are corrected by it)
%   g.pulsesPerScan
%   g.timeS, g.azDeg  time and antenna azimuth of each requested pulse
%                     (azimuth 0 = north, clockwise)

c   = 299792458;
rad = P.radar;
g.cellM   = c / (2 * P.fs * 1e6);
g.rangeM  = rad.rangeOffsetM + (0:nRange-1)' * g.cellM;
g.lambdaM = c / (rad.fcMHz * 1e6);
g.vUnambMps = g.lambdaM * rad.prfHz / 2;

nfft = max(P.fft.nfft, P.fft.nPulses);
k    = (0:nfft-1) - P.fft.shift * floor(nfft/2);       % frequency index per bin
fd   = mod(k/nfft + 0.5, 1) - 0.5;                       % cycles per pulse, -0.5..0.5
g.binVelMps = fd * rad.prfHz * g.lambdaM / 2;

g.degPerPulse   = rad.rpm * 6 / rad.prfHz;
% An output row labelled with pulse p is computed from the pulses before
% it (canceler + FFT window + integration buffer); its centre lags by:
nTaps = numel(rsp_canceler_taps(P.canceler)) * P.canceler.enable + ~P.canceler.enable;
g.outputDelayPulses = ((nTaps - 1) + (P.fft.nPulses - 1) + (P.nci.nFrames - 1) * P.fft.hop) / 2;
g.pulsesPerScan = round(360 / g.degPerPulse);
if nargin >= 3
    pulses  = pulses(:);
    g.timeS = (pulses - 1) / rad.prfHz;
    g.azDeg = mod(rad.azStartDeg + g.degPerPulse * (pulses - 1), 360);
end
end
