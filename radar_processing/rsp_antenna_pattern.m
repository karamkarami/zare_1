function G = rsp_antenna_pattern(dAzDeg, A)
%RSP_ANTENNA_PATTERN One-way azimuth power pattern of the antenna (linear, peak 1).
%
%   G = rsp_antenna_pattern(dAzDeg, A)
%
%   dAzDeg : angle from the beam axis [deg], any size, any range (wrapped)
%   A      : P.antenna
%            .type          'model' (default) | 'sinc' | 'gaussian' | 'table'
%            .beamwidthDeg  3 dB beamwidth                       (4)
%            .sidelobeDb    first sidelobe level, one way          (-25)
%            .backlobeDb    level from 90 deg to the back, one way (-40)
%            .nullDepthDb   depth of the nulls between sidelobes   (-20)
%            .patternDeg    'table': angles from the beam axis [deg], either
%                           -180..180 or 0..180 (symmetric pattern)
%            .patternDb     'table': one-way power gain at those angles [dB]
%                           (normalised to a 0 dB peak)
%            Load a measured pattern from a file with rsp_antenna_load.
%
%   'model'    Main lobe : Gaussian with the given 3 dB width.
%              Sidelobes : lobes one beamwidth apart; the envelope starts at
%              sidelobeDb on the first sidelobe (1.6 beamwidths) and falls
%              (linear in log-angle) to backlobeDb at 90 deg, constant behind.
%   'sinc'     uniformly illuminated aperture, sinc^2 (-13.3 dB sidelobes),
%              never below backlobeDb
%   'gaussian' Gaussian main lobe only, never below backlobeDb
%   'table'    your pattern, interpolated linearly in dB, wrapped in angle;
%              beyond the table its last value
%
%   The two-way echo power of a point target is G^2 (amplitude G).

if ~isfield(A, 'type') || isempty(A.type)
    A.type = 'model';
end
th = mod(dAzDeg + 180, 360) - 180;                    % -180 .. 180 deg
bw = A.beamwidthDeg;
switch lower(A.type)
    case 'model'
        th  = abs(th);
        gm  = exp(-4*log(2) * (th/bw).^2);            % main lobe
        th1 = 1.6 * bw;                               % first sidelobe
        x   = min(max(log(max(th, th1)/th1) / log(90/th1), 0), 1);
        env = A.sidelobeDb + (A.backlobeDb - A.sidelobeDb) * x;      % dB
        rip = cos(pi*(th - th1)/bw).^2;                              % lobes and nulls
        rip = max(rip, 10^(A.nullDepthDb/10));
        gs  = 10.^(env/10) .* rip;
        gs(th < th1 - bw/2) = 0;                      % inside the main lobe
        G   = max(gm, gs);
    case 'sinc'
        x = 0.8859 * th / bw;                         % sinc^2 is 0.5 at x = 0.443
        s = ones(size(x));
        nz = x ~= 0;
        s(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
        G = max(s.^2, 10^(A.backlobeDb/10));
    case 'gaussian'
        G = max(exp(-4*log(2) * (th/bw).^2), 10^(A.backlobeDb/10));
    case 'table'
        a = A.patternDeg(:);
        d = A.patternDb(:);
        if min(a) >= 0                                % 0..180 given: symmetric
            th = abs(th);
        end
        q = min(max(th, min(a)), max(a));             % beyond the table: last value
        G = reshape(10.^(interp1(a, d, q(:), 'linear') / 10), size(th));
    otherwise
        error('rsp_antenna_pattern:type', 'Unknown P.antenna.type "%s".', A.type);
end
end
