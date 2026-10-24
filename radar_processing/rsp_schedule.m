function run = rsp_schedule(P, S, sched, varargin)
%RSP_SCHEDULE Run a schedule of scans: each entry with its own radar settings.
%
%   run = rsp_schedule(P, S, sched, 'option', value, ...)
%
%   P     : base parameters (radar_params)
%   S     : scenario (rsp_simulate); the targets keep moving over the schedule
%   sched : struct array, one entry per scan (or group of revolutions):
%           .name     label                                       ('scan k')
%           .turns    antenna revolutions                         (1)
%           .degrees  rotation instead of turns, e.g. 90 for a quarter
%           .gapS     idle time before the entry [s] (antenna keeps turning,
%                     targets keep moving, nothing is transmitted)   (0)
%           .set      {'path', value, ...} parameter changes for this entry,
%                     e.g. {'radar.prfHz', 1500, 'radar.fcMHz', 1300,
%                           'canceler.order', 2, 'fft.nPulses', 32, 'fft.nfft', 32}
%           .apply    function handle P = f(P) for anything else   ([])
%           Every entry starts from the base P (changes do not carry over).
%
%   Options
%     'display'  true   one PPI for the whole schedule (the sweep continues)
%     'quiet'    false
%     'replay'   false  (unused here; see rsp_ppi_replay on run.scans(k))
%
%   run.scans      result of rsp_scan for every entry (without the display)
%   run.timeline   one row per entry (struct of columns):
%                  entry, name, startS, endS, durationS, startAzDeg, endAzDeg,
%                  turns, prfHz, fcMHz, lambdaM, rangeUnambKm, vUnambMps,
%                  canceler, fftPulses, nciFrames, pulses, found, echoes,
%                  falsePlots, pfa, pfaDesign
%   run.ppi        display state (if displayed)
%
%   Time and azimuth run on from one entry to the next: entry k+1 starts
%   where entry k stopped (+ gapS), so the targets are where they should
%   be. Each entry uses its own noise / clutter seed (S.seed + k - 1).

opt = struct('display', true, 'quiet', false, 'replay', false);
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
if ~isfield(S, 'seed') || isempty(S.seed)
    S.seed = 1;
end
nE   = numel(sched);
tNow = 0;
if isfield(P.radar, 'startTimeS') && ~isempty(P.radar.startTimeS)
    tNow = P.radar.startTimeS;
end
azNow = P.radar.azStartDeg;
ppi   = [];
cols  = {'entry', 'startS', 'endS', 'durationS', 'startAzDeg', 'endAzDeg', 'turns', ...
         'prfHz', 'fcMHz', 'lambdaM', 'rangeUnambKm', 'vUnambMps', 'canceler', ...
         'fftPulses', 'nciFrames', 'pulses', 'found', 'echoes', 'falsePlots', 'pfa', 'pfaDesign'};
for i = 1:numel(cols)
    tl.(cols{i}) = zeros(nE, 1);
end
tl.name = cell(nE, 1);

for k = 1:nE
    e  = sched(k);
    nm = field(e, 'name', sprintf('scan %d', k));
    Q  = P;
    st = field(e, 'set', {});
    for i = 1:2:numel(st)
        Q = setPath(Q, st{i}, st{i+1});
    end
    fn = field(e, 'apply', []);
    if ~isempty(fn)
        Q = fn(Q);
    end
    if Q.fft.nfft < Q.fft.nPulses
        Q.fft.nfft = Q.fft.nPulses;                            % nfft >= nPulses
    end
    deg = field(e, 'degrees', []);
    if isempty(deg)
        deg = 360 * field(e, 'turns', 1);
    end
    gap = field(e, 'gapS', 0);
    degPerS = Q.radar.rpm * 6;
    tNow  = tNow + gap;
    azNow = mod(azNow + degPerS * gap, 360);
    Q.radar.startTimeS = tNow;
    Q.radar.azStartDeg = azNow;
    Sk = S;
    Sk.seed = S.seed + k - 1;

    if ~opt.quiet
        fprintf(['\n===== Schedule %d/%d: %s  (t = %.1f s, az %.1f deg, %.2f turn(s), ' ...
                 'PRF %g Hz, fc %g MHz, %d-pulse canceler, FFT %d, NCI %d) =====\n'], ...
                k, nE, nm, tNow, azNow, deg/360, Q.radar.prfHz, Q.radar.fcMHz, ...
                numel(rsp_canceler_taps(Q.canceler)) * Q.canceler.enable, ...
                Q.fft.nPulses, Q.nci.nFrames);
    end
    sc = rsp_scan(Q, Sk, 'degrees', deg, 'display', opt.display, 'quiet', opt.quiet, ...
                  'ppi', ppi, 'label', sprintf('[%d] %s   ', k, nm));
    if opt.display
        ppi = sc.ppi;
        sc  = rmfield(sc, 'ppi');
    end
    sc.name = nm;
    g  = rsp_geometry(Q, sc.S.nRange);
    nP = numel(sc.azDeg);
    dur = nP / Q.radar.prfHz;

    tl.entry(k)        = k;
    tl.name{k}         = nm;
    tl.startS(k)       = tNow;
    tl.endS(k)         = tNow + dur;
    tl.durationS(k)    = dur;
    tl.startAzDeg(k)   = azNow;
    tl.endAzDeg(k)     = mod(azNow + nP * g.degPerPulse, 360);
    tl.turns(k)        = nP * g.degPerPulse / 360;
    tl.prfHz(k)        = Q.radar.prfHz;
    tl.fcMHz(k)        = Q.radar.fcMHz;
    tl.lambdaM(k)      = g.lambdaM;
    tl.rangeUnambKm(k) = g.rangeUnambM / 1e3;
    tl.vUnambMps(k)    = g.vUnambMps / 2;
    tl.canceler(k)     = numel(rsp_canceler_taps(Q.canceler)) * Q.canceler.enable;
    tl.fftPulses(k)    = Q.fft.nPulses;
    tl.nciFrames(k)    = Q.nci.nFrames;
    tl.pulses(k)       = nP;
    tl.found(k)        = sc.falseAlarms.found;
    tl.echoes(k)       = sc.falseAlarms.echoes;
    tl.falsePlots(k)   = sc.falseAlarms.falsePlots(end);
    tl.pfa(k)          = sc.falseAlarms.pfa(end);
    tl.pfaDesign(k)    = sc.falseAlarms.pfaDesign;
    scans(k)           = sc;                                    %#ok<AGROW>

    tNow  = tNow + dur;
    azNow = tl.endAzDeg(k);
end

run.scans    = scans;
run.timeline = tl;
run.P        = P;
run.S        = S;
run.sched    = sched;
if opt.display
    run.ppi = ppi;
end
if ~opt.quiet
    rsp_schedule_print(tl);
end
end

% ======================================================================================
function v = field(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function P = setPath(P, path, value)
% P.a.b = value for path 'a.b'; indices like 'pulse(2).fcMHz' are allowed.
parts = strsplit(path, '.');
sub   = struct('type', {}, 'subs', {});
for i = 1:numel(parts)
    t = regexp(parts{i}, '^(\w+)\((\d+)\)$', 'tokens', 'once');
    if isempty(t)
        sub(end+1) = struct('type', '.', 'subs', parts{i});          %#ok<AGROW>
    else
        sub(end+1) = struct('type', '.', 'subs', t{1});              %#ok<AGROW>
        sub(end+1) = struct('type', '()', 'subs', {{str2double(t{2})}});  %#ok<AGROW>
    end
end
P = subsasgn(P, sub, value);
end
