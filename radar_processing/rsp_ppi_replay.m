function h = rsp_ppi_replay(scan, varargin)
%RSP_PPI_REPLAY Replay a processed scan on the PPI with a rotating sweep.
%
%   h = rsp_ppi_replay(scan, 'option', value, ...)
%
%   scan : result of rsp_scan
%
%   Options
%     'speed'    1      1 = real antenna speed (P.radar.rpm), 4 = four times faster,
%                       Inf = as fast as possible
%     'stepDeg'  1      sweep step per screen update [deg]
%     'bins'     []     Doppler bins shown ([] = as processed, P.output.bins);
%                       a subset of the processed bins, e.g. 5 or [2:16]
%     'P'        []     parameters for the display (default scan.P)
%
%   The output after the CFAR (max over bins) is painted; plots are shown
%   when the sweep passes them.

opt = struct('speed', 1, 'stepDeg', 1, 'bins', [], 'P', []);
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
P = scan.P;
if ~isempty(opt.P)
    P = opt.P;
end
P.ppi.source = 'max';
g     = rsp_geometry(P, numel(scan.rangeM));
nPul  = numel(scan.azDeg);
step  = max(1, round(opt.stepDeg / g.degPerPulse));

V = scan.maxValue;
plots = scan.plots;
if ~isempty(opt.bins)                                   % keep only the chosen bins
    keep = ismember(scan.det.bin, opt.bins);
    V = sparse(scan.det.pulse(keep), scan.det.cell(keep), scan.det.value(keep), ...
               nPul, numel(scan.rangeM));
    d = subset(scan.det, keep);
    plots = rsp_extract_plots(d, P, numel(scan.rangeM));
end
plotPulse = round(plots.pulse + g.outputDelayPulses);  % when the sweep reaches them

h = rsp_ppi_init(P, scan.rangeM, 'PPI replay');
t0 = tic;
for p0 = 1:step:nPul
    p    = (p0:min(p0 + step - 1, nPul))';
    rows = 10*log10(full(V(p, :)));
    sel  = plotPulse >= p(1) & plotPulse <= p(end);
    h = rsp_ppi_update(h, scan.azDeg(p), rows, subset(plots, sel), ...
        sprintf('Az %5.1f deg   t = %5.1f s   bins %s', scan.azDeg(p(end)), ...
                (p(end) - 1) / P.radar.prfHz, binText(opt.bins, P.output.bins)));
    if isfinite(opt.speed)
        wait = (p(end) / P.radar.prfHz) / opt.speed - toc(t0);
        if wait > 0
            pause(wait);
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
