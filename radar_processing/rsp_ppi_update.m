function h = rsp_ppi_update(h, azDeg, rowsDb, plots, info)
%RSP_PPI_UPDATE Paint newly received pulses on the PPI and move the sweep.
%
%   h = rsp_ppi_update(h, azDeg, rowsDb, plots, info)
%
%   h      : state from rsp_ppi_init
%   azDeg  : antenna azimuth of each new row [deg] (increasing, clockwise)
%   rowsDb : rows x range cells, values in dB (-Inf = nothing)
%   plots  : new plots to mark (struct with azDeg, rangeM, velocityMps), or []
%   info   : text for the title (optional)
%
%   Range is reduced to the display resolution with a peak hold (a
%   detection narrower than a pixel stays visible); this is display only.

k = numel(azDeg);
if k == 0
    return
end
azDeg = azDeg(:);

% --- unwrapped azimuth of the new rows ---------------------------------------------------
if isempty(h.azAbs)
    h.azAbs = azDeg(1) - 1e-6;
end
step = mod(diff([h.azAbs; azDeg]) + 180, 360) - 180;
a    = h.azAbs + cumsum(step);
a0   = h.azAbs;
span = min(a(end) - a0, 360);

% --- rows at display resolution (peak hold over range) --------------------------------
nc = min(size(rowsDb, 2), h.nCells);
r  = -inf(k, h.nGrp * h.dec);
r(:, 1:nc) = rowsDb(:, 1:nc);
r  = squeeze(max(reshape(r', h.dec, h.nGrp, k), [], 1))';
if k == 1
    r = r(:)';
end

if isempty(h.clim)
    v = r(isfinite(r));
    if isempty(v)
        h.clim = [0 30];
    elseif strcmpi(h.P.ppi.source, 'max')
        h.clim = [min(v) max(min(v) + 25, max(v))];
    else
        h.clim = rsp_db_limits(v, 60, 99.99);
    end
end

% --- paint the swept sector --------------------------------------------------------------------
dd  = mod(h.azSorted - a0, 360);
sel = dd > 0 & dd <= span;
if span >= 360
    sel(:) = true;
end
if any(sel)
    [au, iu] = unique(a - a0);
    row = interp1(au, iu, dd(sel), 'nearest', 'extrap');
    h.img(h.pix(sel))   = r(sub2ind(size(r), row, h.grp(sel)));
    h.paint(h.pix(sel)) = a0 + dd(sel);
end
h.azAbs = a(end);

% --- afterglow and screen update --------------------------------------------------------------
age  = h.azAbs - h.paint;                                 % 0 .. 360 deg
img  = h.img - h.fade * age / 360;
img(~isfinite(img)) = h.clim(1) - 1;
set(h.im, 'CData', img);
set(h.ax, 'CLim', h.clim);
th = mod(h.azAbs, 360) * pi/180;
set(h.sweep, 'XData', [0 sin(th)] * h.rMax/1e3, 'YData', [0 cos(th)] * h.rMax/1e3);

% --- plots: replace the ones swept again, add the new ones -----------------------------------
if h.P.ppi.showPlots
    old = (h.azAbs - h.plots.az) >= 360;
    h.plots.az(old) = [];  h.plots.x(old) = [];  h.plots.y(old) = [];  h.plots.v(old) = [];
    if nargin >= 4 && ~isempty(plots) && ~isempty(plots.azDeg)
        pa = h.azAbs - mod(h.azAbs - plots.azDeg(:), 360);          % unwrapped
        h.plots.az = [h.plots.az; pa];
        h.plots.x  = [h.plots.x;  plots.rangeM(:) .* sin(plots.azDeg(:)*pi/180) / 1e3];
        h.plots.y  = [h.plots.y;  plots.rangeM(:) .* cos(plots.azDeg(:)*pi/180) / 1e3];
        h.plots.v  = [h.plots.v;  plots.velocityMps(:)];
    end
    set(h.marker, 'XData', h.plots.x, 'YData', h.plots.y);
    if ~isempty(h.labels)
        delete(h.labels(ishandle(h.labels)));
    end
    h.labels = zeros(numel(h.plots.x), 1);
    for i = 1:numel(h.plots.x)
        h.labels(i) = text(h.plots.x(i) + 0.012*h.rMax/1e3, h.plots.y(i), ...
            sprintf('%.0f m/s', h.plots.v(i)), 'Color', [1 0.85 0.2], ...
            'FontSize', 8, 'Parent', h.ax);
    end
end

if nargin >= 5 && ~isempty(info)
    set(h.title, 'String', info);
end
drawnow;
end
