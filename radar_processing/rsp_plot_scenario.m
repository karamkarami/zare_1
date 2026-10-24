function fig = rsp_plot_scenario(P, S, varargin)
%RSP_PLOT_SCENARIO Scenario at a glance: clutter map, targets and antenna pattern.
%
%   fig = rsp_plot_scenario(P, S)
%   fig = rsp_plot_scenario(P, S, 'visible', 'off')
%
%   Left : mean clutter power of every layer (dB, same scale as S.noiseDb),
%          terrain shadow included, the region outlines and the targets
%          (number, range, velocity) at their starting position
%   Right: one-way and two-way antenna pattern, full circle and main lobe,
%          with the 3 dB beamwidth

opt = struct('visible', 'on');
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
P.antenna = rsp_antenna_load(P.antenna);
S  = rsp_levels(S);
g  = rsp_geometry(P, S.nRange);
L  = rsp_clutter_types([]);
if isfield(S, 'clutter')
    L = rsp_clutter_types(S.clutter);
end
rMax = g.rangeM(end);

fig = figure('Name', 'Scenario', 'Color', 'w', 'Position', [60 60 1300 620], ...
             'Visible', opt.visible);

% --- clutter map --------------------------------------------------------------------
ax = subplot(1, 2, 1);
azG = (0:0.5:359.5)';
rG  = (0:250:rMax);
pw  = zeros(numel(azG), numel(rG));
for k = 1:numel(L)
    pw = pw + 10^(L(k).powerDb/10) * rsp_clutter_mask(L(k), azG, rG) .* ...
         rsp_clutter_shadow(L, azG, rG, k);
end
db  = 10*log10(pw + 10^(S.noiseDb/10));                % noise floor where no clutter
N   = 500;
x   = linspace(-rMax, rMax, N);
[X, Y] = meshgrid(x, x);
ia  = mod(round(mod(atan2d(X, Y), 360) / 0.5), numel(azG)) + 1;
ir  = round(hypot(X, Y) / 250) + 1;
img = S.noiseDb * ones(N);                             % outside the record: noise level
in  = ir <= numel(rG);
img(in) = db(sub2ind(size(db), ia(in), ir(in)));
imagesc(x/1e3, x/1e3, img, 'Parent', ax);
set(ax, 'YDir', 'normal');
axis(ax, 'equal');
axis(ax, [-1 1 -1 1] * rMax/1e3);
colormap(ax, flipud(gray(256)) * 0.85 + 0.15);
cb = colorbar(ax);
ylabel(cb, 'mean clutter power per sample [dB]');
hold(ax, 'on');
col = struct('ground', [0.6 0.4 0.2], 'mountain', [0.35 0.35 0.45], 'forest', [0.1 0.55 0.1]);
for k = 1:numel(L)
    [ox, oy] = rsp_clutter_outline(L(k));
    plot(ax, ox, oy, '-', 'Color', col.(L(k).type), 'LineWidth', 1.5);
    [~, j] = max(ox.^2 + oy.^2);
    text(ox(j), oy(j), sprintf(' %s %.0f dB', L(k).name, L(k).powerDb), 'Parent', ax, ...
         'Color', col.(L(k).type), 'FontSize', 9, 'FontWeight', 'bold');
end
rings(ax, rMax, P.ppi.ringKm);
for k = 1:numel(S.targets)
    t = S.targets(k);
    tx = t.rangeM * sind(t.azDeg) / 1e3;
    ty = t.rangeM * cosd(t.azDeg) / 1e3;
    plot(ax, tx, ty, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    text(tx, ty, sprintf('  T%d  %.0f km  %.0f m/s  %.0f dB', k, t.rangeM/1e3, ...
         t.velocityMps, t.powerDb), 'Parent', ax, 'Color', [0.7 0 0], 'FontSize', 8);
end
xlabel(ax, 'East [km]');  ylabel(ax, 'North [km]');
title(ax, sprintf(['Scenario: clutter and targets (noise %.0f dB, PRF %g Hz, %g MHz, ' ...
      '%g rpm)'], S.noiseDb, P.radar.prfHz, P.radar.fcMHz, P.radar.rpm));

% --- antenna pattern ----------------------------------------------------------------
th = (-180:0.05:180)';
G  = 10*log10(rsp_antenna_pattern(th, P.antenna));
bw = P.antenna.beamwidthDeg;
ax = subplot(2, 2, 2);
plot(ax, th, G, 'b-', th, 2*G, 'r-');
grid(ax, 'on');
axis(ax, [-180 180 max(min(2*G), -100) 3]);
set(ax, 'XTick', -180:45:180);
legend(ax, 'one way', 'two way (echo)', 'Location', 'south');
ylabel(ax, 'gain [dB]');
title(ax, sprintf('Antenna pattern (%s): 3 dB beamwidth %.2f deg', P.antenna.type, bw));
ax = subplot(2, 2, 4);
k  = abs(th) <= 5*bw;
plot(ax, th(k), G(k), 'b-', th(k), 2*G(k), 'r-', [-bw bw]/2, [-3 -3], 'ko-');
grid(ax, 'on');
axis(ax, [-5*bw 5*bw -60 3]);
xlabel(ax, 'angle from the beam axis [deg]');
ylabel(ax, 'gain [dB]');
title(ax, 'Main lobe and first sidelobes');
end

function rings(ax, rMax, ringKm)
t = linspace(0, 2*pi, 361);
for r = ringKm:ringKm:rMax/1e3
    plot(ax, r*sin(t), r*cos(t), ':', 'Color', [0.5 0.5 0.5]);
end
end
