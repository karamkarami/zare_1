function h = rsp_ppi_replay(scan, varargin)
%RSP_PPI_REPLAY Replay a processed scan on the PPI with a rotating sweep.
%
%   h = rsp_ppi_replay(scan, 'option', value, ...)
%
%   scan : result of rsp_scan
%
%   Options
%     'speed'    P.ppi.speed    1 = real antenna speed (P.radar.rpm), 4 = four
%                               times faster, Inf = as fast as possible
%     'stepDeg'  P.ppi.stepDeg  sweep step per screen update [deg]
%     'loops'    1              how many times the recording is played
%     'bins'     []             Doppler bins shown ([] = as processed,
%                               P.output.bins), e.g. 5 or [2:16]
%     'P'        []             parameters for the display (default scan.P)
%     'framesDir'  ''           folder: save a PNG of the screen at every
%                               'frameEvery'-th update (to make a video)
%     'frameEvery' 1
%
%   The output after the CFAR is painted (as dB over the threshold, or as
%   the value: P.ppi.maxScale); each plot appears when the sweep has passed
%   it, the plots of earlier revolutions stay as a trail.

P = scan.P;
opt = struct('speed', P.ppi.speed, 'stepDeg', P.ppi.stepDeg, 'loops', 1, 'bins', [], 'P', [], ...
             'framesDir', '', 'frameEvery', 1);
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
if ~isempty(opt.P)
    P = opt.P;
end
P.ppi.source = 'max';
R     = numel(scan.rangeM);
g     = rsp_geometry(P, R);
nPul  = numel(scan.azDeg);
step  = max(1, round(opt.stepDeg / g.degPerPulse));
lag   = ceil(1.5 * P.antenna.beamwidthDeg / g.degPerPulse);

det = scan.det;
if ~isempty(opt.bins)                                   % keep only the chosen bins
    keep = ismember(det.bin, opt.bins);
    f = fieldnames(det);
    for i = 1:numel(f)
        det.(f{i}) = det.(f{i})(keep);
    end
end
if strcmpi(P.ppi.maxScale, 'margin')
    V = sparse(det.pulse, det.cell, max(det.marginDb, 1e-3), nPul, R);
else
    V = sparse(det.pulse, det.cell, 10*log10(det.value), nPul, R);
end
plots = rsp_extract_plots(det, P, R);
showAt = plots.pulse + lag;                             % shown once the beam has passed

h = rsp_ppi_init(P, scan.rangeM, 'PPI replay');
nUpd = 0;
if ~isempty(opt.framesDir) && ~exist(opt.framesDir, 'dir')
    mkdir(opt.framesDir);
end
for loop = 1:opt.loops
    t0 = tic;
    for p0 = 1:step:nPul
        p    = (p0:min(p0 + step - 1, nPul))';
        rows = full(V(p, :));
        rows(rows == 0) = -Inf;
        sel  = showAt >= p(1) & showAt <= p(end) | (p(end) == nPul & showAt > nPul);
        h = rsp_ppi_update(h, scan.azDeg(p), rows, subset(plots, sel), ...
            sprintf('Az %5.1f deg   t = %5.1f s   scan %d   bins %s', scan.azDeg(p(end)), ...
                    (p(end) - 1) / P.radar.prfHz, floor((p(end) - 1) * g.degPerPulse / 360) + 1, ...
                    binText(opt.bins, P.output.bins)));
        nUpd = nUpd + 1;
        if ~isempty(opt.framesDir) && mod(nUpd - 1, opt.frameEvery) == 0
            print(h.fig, fullfile(opt.framesDir, sprintf('ppi_%05d.png', nUpd)), '-dpng', '-r72');
        end
        if isfinite(opt.speed)
            wait = (p(end) / P.radar.prfHz) / opt.speed - toc(t0);
            if wait > 0
                pause(wait);
            end
        end
    end
end
end

function s = subset(s, k)
f = fieldnames(s);
for i = 1:numel(f)
    s.(f{i}) = s.(f{i})(k);
end
end

function t = binText(b, b0)
if isempty(b)
    b = b0;
end
if isempty(b)
    t = 'all';
else
    t = mat2str(b);
end
end
