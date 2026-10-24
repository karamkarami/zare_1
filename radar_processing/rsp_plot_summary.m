function fig = rsp_plot_summary(scan, varargin)
%RSP_PLOT_SUMMARY One scan at a glance: detections, plots, truth, false alarms.
%
%   fig = rsp_plot_summary(scan)
%   fig = rsp_plot_summary(scan, 'visible', 'off', 'maxPoints', 30000)
%
%   Left : every output detection (coloured by its margin over the CFAR
%          threshold), the clutter regions, the true target crossings
%          (circles), the plots of targets (green squares) and the false
%          plots (red crosses)
%   Right: measured false alarm probability per clutter region against the
%          design value, and the target table (found, errors, margin)

opt = struct('visible', 'on', 'maxPoints', 30000);
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
P  = scan.P;
fa = scan.falseAlarms;
g  = rsp_geometry(P, numel(scan.rangeM));
rMax = scan.rangeM(end) / 1e3;
nm = '';
if isfield(scan, 'name')
    nm = [scan.name ': '];
end

fig = figure('Name', 'Scan summary', 'Color', 'w', 'Position', [60 60 1300 620], ...
             'Visible', opt.visible);

% --- map ----------------------------------------------------------------------------
ax = subplot(1, 2, 1);
hold(ax, 'on');
d  = scan.det;
az = mod(d.azDeg - g.degPerPulse * g.outputDelayPulses, 360);
k  = 1:numel(az);
if numel(k) > opt.maxPoints
    k = round(linspace(1, numel(az), opt.maxPoints));
end
if ~isempty(k)
    scatter(ax, d.rangeM(k) .* sind(az(k)) / 1e3, d.rangeM(k) .* cosd(az(k)) / 1e3, 4, ...
            min(d.marginDb(k), 40), 'filled');
end
colormap(ax, jet(256));
cb = colorbar(ax);
ylabel(cb, 'margin over the CFAR threshold [dB]');
col = struct('ground', [0.6 0.4 0.2], 'mountain', [0.35 0.35 0.45], 'forest', [0.1 0.55 0.1]);
for i = 1:numel(scan.clutter)
    [ox, oy] = rsp_clutter_outline(scan.clutter(i));
    plot(ax, ox, oy, '-', 'Color', col.(scan.clutter(i).type), 'LineWidth', 1.5);
end
t = linspace(0, 2*pi, 361);
for r = P.ppi.ringKm:P.ppi.ringKm:rMax
    plot(ax, r*sin(t), r*cos(t), ':', 'Color', [0.6 0.6 0.6]);
end
tr = scan.truth;
if ~isempty(tr)
    plot(ax, [tr.apparentRangeM] .* sind([tr.azDeg]) / 1e3, ...
         [tr.apparentRangeM] .* cosd([tr.azDeg]) / 1e3, 'ko', 'MarkerSize', 12, 'LineWidth', 1.2);
end
rep = scan.report;
fnd = rep([rep.found]);
if ~isempty(fnd)
    plot(ax, [fnd.plotRangeM] .* sind([fnd.plotAzDeg]) / 1e3, ...
         [fnd.plotRangeM] .* cosd([fnd.plotAzDeg]) / 1e3, 'gs', 'MarkerSize', 9, ...
         'LineWidth', 1.5);
end
if ~isempty(fa.plots.azDeg)
    plot(ax, fa.plots.rangeM .* sind(fa.plots.azDeg) / 1e3, ...
         fa.plots.rangeM .* cosd(fa.plots.azDeg) / 1e3, 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
end
axis(ax, 'equal');
axis(ax, [-1 1 -1 1] * rMax);
box(ax, 'on');
xlabel(ax, 'East [km]');  ylabel(ax, 'North [km]');
title(ax, {sprintf('%sdetections, target plots (green), truth (o), false plots (x)', nm), ...
      sprintf('PRF %g Hz, %g MHz, %d-pulse canceler, FFT %d, NCI %d, bins %s', P.radar.prfHz, ...
      P.radar.fcMHz, numel(rsp_canceler_taps(P.canceler)) * P.canceler.enable, ...
      P.fft.nPulses, P.nci.nFrames, binText(P.output.bins))}, 'Interpreter', 'none', ...
      'FontSize', 9);

% --- false alarm probability per region ---------------------------------------------
ax = subplot(2, 2, 2);
sel = find(fa.cells > 0);
v   = log10(fa.pfa(sel));
v(~isfinite(v)) = NaN;                                  % no false detection
bar(ax, v);
hold(ax, 'on');
for i = find(isnan(v(:)'))
    text(i, -1, '0', 'Parent', ax, 'HorizontalAlignment', 'center');
end
if isfinite(fa.pfaDesign)
    plot(ax, [0.5 numel(sel) + 0.5], log10(fa.pfaDesign) * [1 1], 'r--', 'LineWidth', 1.5);
end
ylim(ax, [floor(min([v(:); log10(fa.pfaDesign); -6])) - 1, 0]);
set(ax, 'XTick', 1:numel(sel), 'XTickLabel', fa.names(sel));
ylabel(ax, 'log_{10} Pfa per output cell');
title(ax, sprintf('False alarms: measured (bars), design (red): %d false plots, %.1f per scan', ...
      fa.falsePlots(end), fa.plotsPerScan(end)));
grid(ax, 'on');

% --- target table -------------------------------------------------------------------
ax = subplot(2, 2, 4);
axis(ax, 'off');
lines = {sprintf('%-4s %-4s %8s %8s %8s  %-9s %7s %7s %7s %7s', 'tgt', 'scan', 'az', ...
         'range', 'vel', 'result', 'dAz', 'dR', 'dV', 'margin')};
for i = 1:min(numel(rep), 14)
    r = rep(i);
    if r.found
        lines{end+1} = sprintf('%-4d %-4d %8.1f %8.0f %8.1f  %-9s %7.2f %7.0f %7.1f %7.1f', ...
            r.target, r.scan, r.azDeg, r.rangeM, r.velocityMps, 'found', r.dAzDeg, ...
            r.dRangeM, r.dVelMps, r.marginDb);                          %#ok<AGROW>
    else
        why = 'MISSED';
        if r.shadowDb < -3, why = 'SHADOWED'; end
        lines{end+1} = sprintf('%-4d %-4d %8.1f %8.0f %8.1f  %-9s', r.target, r.scan, ...
            r.azDeg, r.rangeM, r.velocityMps, why);                     %#ok<AGROW>
    end
end
lines{end+1} = sprintf('Pd %.2f (%d of %d echoes), measured Pfa %.1e (design %.1e)', fa.pd, ...
                       fa.found, fa.echoes, fa.pfa(end), fa.pfaDesign);
text(0, 1, lines, 'Parent', ax, 'FontName', 'FixedWidth', 'FontSize', 8, ...
     'VerticalAlignment', 'top', 'Interpreter', 'none');
end

function t = binText(b)
if isempty(b)
    t = 'all';
else
    t = mat2str(b);
end
end
