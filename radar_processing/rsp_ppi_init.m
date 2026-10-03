function h = rsp_ppi_init(P, rangeM, ttl)
%RSP_PPI_INIT Create a PPI (plan position indicator) display.
%
%   h = rsp_ppi_init(P, rangeM, ttl)
%
%   P      : parameters (P.ppi, P.radar)
%   rangeM : range of every cell [m] (rsp_geometry(...).rangeM)
%   ttl    : figure title (optional)
%   h      : display state, pass it to rsp_ppi_update
%
%   North up, azimuth clockwise. Each update paints the sector swept since
%   the previous update; older echoes fade by P.ppi.fadeDb per revolution
%   (afterglow), like a phosphor screen. Range is shown in km.

if nargin < 3
    ttl = 'PPI';
end
Q      = P.ppi;
rangeM = rangeM(:);
cellM  = rangeM(2) - rangeM(1);
if isempty(Q.maxRangeKm)
    rMax = rangeM(end);
else
    rMax = min(Q.maxRangeKm * 1e3, rangeM(end));
end
N = Q.pixels;

% --- pixel geometry (computed once) --------------------------------------------------
x       = linspace(-rMax, rMax, N);
[X, Y]  = meshgrid(x, x);
rho     = hypot(X, Y);
az      = mod(atan2(X, Y) * 180/pi, 360);                % 0 = north, clockwise
inside  = rho <= rMax & rho >= rangeM(1);
nCells  = find(rangeM <= rMax, 1, 'last');
dec     = max(1, floor(nCells / (N/2)));                 % range cells per display bin
cellPix = min(nCells, floor((rho(inside) - rangeM(1)) / cellM) + 1);
[h.azSorted, ord] = sort(az(inside));
pix     = find(inside);
h.pix   = pix(ord);
h.grp   = ceil(cellPix(ord) / dec);
h.dec   = dec;
h.nCells = nCells;
h.nGrp  = ceil(nCells / dec);
h.N     = N;
h.img   = nan(N, N);                                     % painted value [dB]
h.paint = nan(N, N);                                     % unwrapped azimuth of the paint
h.azAbs = [];                                            % unwrapped sweep azimuth
h.clim  = Q.climDb;
h.fade  = Q.fadeDb;
h.P     = P;
h.plots = struct('az', zeros(0, 1), 'x', zeros(0, 1), 'y', zeros(0, 1), 'v', zeros(0, 1));

% --- figure -------------------------------------------------------------------------------
h.fig = figure('Name', ttl, 'Color', 'k', 'Position', [80 60 900 900]);
h.ax  = axes('Parent', h.fig, 'Color', 'k', 'Position', [0.06 0.06 0.84 0.86]);
xk    = x / 1e3;
h.im  = imagesc(xk, xk, zeros(N), 'Parent', h.ax);
set(h.ax, 'YDir', 'normal', 'XColor', [0 0.6 0], 'YColor', [0 0.6 0]);
axis(h.ax, 'equal');
axis(h.ax, [-1 1 -1 1] * rMax/1e3);
hold(h.ax, 'on');
if strcmpi(Q.colormap, 'phosphor')
    colormap(h.ax, phosphor(256));
else
    colormap(h.ax, Q.colormap);
end
cb = colorbar(h.ax);
try
    set(cb, 'Color', [0 0.8 0]);                         % MATLAB colorbar object
catch
    set(cb, 'XColor', [0 0.8 0], 'YColor', [0 0.8 0]);   % Octave colorbar axes
end
ylabel(cb, 'dB');
xlabel(h.ax, 'East [km]');
ylabel(h.ax, 'North [km]');

% range rings and azimuth spokes
gr = [0 0.45 0];
t  = linspace(0, 2*pi, 361);
for r = Q.ringKm:Q.ringKm:rMax/1e3
    plot(h.ax, r*sin(t), r*cos(t), ':', 'Color', gr);
    text(r*sin(pi/36), r*cos(pi/36), sprintf('%g', r), 'Color', gr, 'Parent', h.ax, ...
         'FontSize', 8);
end
for a = 0:30:330
    plot(h.ax, [0 sin(a*pi/180)] * rMax/1e3, [0 cos(a*pi/180)] * rMax/1e3, ':', 'Color', gr);
    text(1.04*rMax/1e3*sin(a*pi/180), 1.04*rMax/1e3*cos(a*pi/180), sprintf('%d', a), ...
         'Color', [0 0.8 0], 'HorizontalAlignment', 'center', 'Parent', h.ax, 'FontSize', 9);
end
h.sweep  = plot(h.ax, [0 0], [0 rMax/1e3], '-', 'Color', [0.6 1 0.6], 'LineWidth', 1.5);
h.trail  = plot(h.ax, NaN, NaN, '.', 'Color', [0.75 0.6 0.15], 'MarkerSize', 8);
h.marker = plot(h.ax, NaN, NaN, 's', 'Color', [1 0.85 0.2], 'MarkerSize', 10, 'LineWidth', 1.2);
h.labels = {};
h.labelSig = [];
h.title  = title(h.ax, ttl, 'Color', [0.6 1 0.6]);
h.rMax   = rMax;
end

function c = phosphor(n)
% Black -> green -> yellow-white, like a radar phosphor screen.
t = linspace(0, 1, n)';
c = [max(0, (t - 0.65) / 0.35).^1.5, min(1, t / 0.55).^0.9, max(0, (t - 0.8) / 0.2) * 0.8];
c(1, :) = 0;
c = min(max(c, 0), 1);
end
