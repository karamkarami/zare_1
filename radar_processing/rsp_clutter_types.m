function L = rsp_clutter_types(clutter)
%RSP_CLUTTER_TYPES Clutter layers of a scenario, with the defaults of each terrain type.
%
%   L = rsp_clutter_types(S.clutter)
%   T = rsp_clutter_types('ground')     defaults of one type (also 'mountain', 'forest')
%
%   S.clutter is a struct array, one layer per clutter region. Every field
%   except .type is optional and takes the default of its type:
%
%     type        'ground' | 'mountain' | 'forest'
%     name        label in reports (default: the type)
%     powerDb     mean clutter power per sample at the receiver [dB] (same
%                 scale as S.noiseDb); for 'mountain' the power at the ridge
%     azDeg       [first last] azimuth sector, clockwise [deg]; [0 360] = all
%                 (may wrap: [350 20] spans 30 deg through north)
%     rangeM      [near far] range extent [m]
%     edgeDeg     soft edge in azimuth [deg]        (power ramps up inside)
%     edgeM       soft edge in range [m]
%     sigmaVMps   spread of the internal motion (wind on trees, ...) [m/s]
%     dcFraction  share of the power without internal motion (rocks, trunks)
%     textureDb   spatial fluctuation of the mean power (log-normal) [dB]
%     slopeDb     power rise from the near to the far edge [dB] (a slope
%                 facing the radar: grazing angle grows towards the ridge)
%     shadow      true = hides what lies behind it (larger range, same azimuth)
%     shadowDb    two-way loss of the echoes in its shadow [dB]
%     scanSpread  true = add the spectral spread of the antenna scanning
%                 (0.265 * PRF / pulses per beamwidth) to every component
%
%   Defaults (noise at 0 dB):
%                 powerDb  sigmaVMps  dcFraction  textureDb  slopeDb  shadow
%     ground         10       0.1        0.9         3          0       no
%     mountain       45       0.05       0.98        6         10       yes (-40 dB)
%     forest         20       0.6        0.5         4          0       no
%
%   The single-struct form of the earlier versions (powerDb, maxRangeM,
%   sigmaVMps, textureDb, without .type) is one 'ground' layer over all
%   azimuths, from 0 to maxRangeM, fading over the last 30 %, with exactly
%   the given spread (no dc part, no scanning spread).

if ischar(clutter)
    L = defaults(clutter);
    return
end
L = struct('type', {}, 'name', {}, 'powerDb', {}, 'azDeg', {}, 'rangeM', {}, ...
           'edgeDeg', {}, 'edgeM', {}, 'sigmaVMps', {}, 'dcFraction', {}, ...
           'textureDb', {}, 'slopeDb', {}, 'shadow', {}, 'shadowDb', {}, 'scanSpread', {});
if isempty(clutter)
    return
end
f = fieldnames(L);
for k = 1:numel(clutter)
    c = clutter(k);
    if ~isfield(c, 'type') || isempty(c.type)
        % earlier single-layer form
        d = defaults('ground');
        d.name       = 'ground';
        d.azDeg      = [0 360];
        d.rangeM     = [0 c.maxRangeM];
        d.edgeM      = 0.3 * c.maxRangeM;
        d.dcFraction = 0;
        d.scanSpread = false;
        d.slopeDb    = 0;
        d.shadow     = false;
        g = {'powerDb', 'sigmaVMps', 'textureDb'};
        for i = 1:numel(g)
            if isfield(c, g{i}) && ~isempty(c.(g{i}))
                d.(g{i}) = c.(g{i});
            end
        end
    else
        d = defaults(c.type);
        g = fieldnames(c);
        for i = 1:numel(g)
            if isfield(d, g{i}) && ~isempty(c.(g{i}))
                d.(g{i}) = c.(g{i});
            end
        end
    end
    if isempty(d.name)
        d.name = d.type;
    end
    for i = 1:numel(f)
        L(k).(f{i}) = d.(f{i});
    end
end
end

% ======================================================================================
function d = defaults(type)
d = struct('type', lower(type), 'name', '', 'powerDb', 10, 'azDeg', [0 360], ...
           'rangeM', [0 25000], 'edgeDeg', 1, 'edgeM', 500, 'sigmaVMps', 0.1, ...
           'dcFraction', 0.9, 'textureDb', 3, 'slopeDb', 0, 'shadow', false, ...
           'shadowDb', -40, 'scanSpread', true);
switch lower(type)
    case 'ground'
    case 'mountain'
        d.powerDb    = 45;
        d.azDeg      = [100 140];
        d.rangeM     = [40000 55000];
        d.sigmaVMps  = 0.05;
        d.dcFraction = 0.98;
        d.textureDb  = 6;
        d.slopeDb    = 10;
        d.shadow     = true;
        d.edgeM      = 1000;
    case 'forest'
        d.powerDb    = 20;
        d.azDeg      = [200 250];
        d.rangeM     = [12000 30000];
        d.sigmaVMps  = 0.6;
        d.dcFraction = 0.5;
        d.textureDb  = 4;
    otherwise
        error('rsp_clutter_types:type', ...
              'Unknown clutter type "%s" (ground, mountain, forest).', type);
end
end
