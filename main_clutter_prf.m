% Clutter-map MTI with a scan-to-scan PRF set, on simulated data.
%
%   clutter_prf : the processor (does not know where the data comes from)
%   clutter_sim : the test-data simulator (does not know the processor)
%
% The radar is stationary for the first scans, then starts to move. Each scan
% prints the PRF, the clutter before / after the filter, and every target's
% power after the filter next to the power expected from the gain of the
% filter used in that cell at the target's velocity (input SNR + beam loss +
% range straddle loss + gain).

if ~exist('doPlot', 'var'), doPlot = true; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'clutter_prf'), fullfile(here, 'clutter_sim'));

%% Processor settings (your radar)
cfg = cpf_config();
cfg.radar.fcHz         = 600e6;
cfg.radar.beamwidthDeg = 4;
cfg.radar.rpm          = 6;
cfg.prf.minHz          = 800;
cfg.prf.maxHz          = 1500;
cfg.prf.count          = 4;        % P
cfg.mti.maxLength      = 8;        % maximum filter length
cfg.target.vMinMps     = 20;
cfg.target.vMaxMps     = 350;

%% Simulator settings (the same radar, written separately on purpose)
S = csim_config();
S.fcHz = cfg.radar.fcHz;  S.beamwidthDeg = cfg.radar.beamwidthDeg;  S.rpm = cfg.radar.rpm;
S.swerling = false;                % steady targets, so the powers can be compared with the design
nScans    = 12;
moveScan  = 9;                     % the platform moves from this scan on
moveVelEN = [10.6 10.6];           % 15 m/s towards north-east

%% Run
st  = cpf_init(cfg);               % first PRF set, from cfg.clutter
cpf_report(st.sel);
cfgNoComp = cfg;  cfgNoComp.mti.motionComp = false;
stNoComp  = [];                    % same processor without motion compensation (moving part)
sim = csim_init(S);
lam = 299792458 / cfg.radar.fcHz;
nT  = numel(S.targets);
th  = linspace(-S.beamwidthDeg / 4, S.beamwidthDeg / 4, 101);
beamLossDb = 10 * log10(mean(exp(-8 * log(2) * th.^2 / S.beamwidthDeg^2)));  % pulses within +-bw/4
resid = nan(nScans, 3);            % clutter in, out, out without motion compensation [dB re noise]
tgtDb = nan(nScans, nT);
tgtPred = nan(nScans, nT);
prfUsed = nan(nScans, 1);

fprintf('\nscan  PRF   platform   clutter in  out (no comp)   target power after the filter [dB re noise] (expected)\n');
for scan = 1:nScans
    if scan == moveScan
        sim.vel = moveVelEN;
        stNoComp = st;  stNoComp.cfg = cfgNoComp;
    end
    [prf, st, sel] = cpf_next_prf(st, norm(sim.vel));
    if ~isempty(sel)
        if sel.changed, fprintf('--- new PRF set\n'); else, fprintf('--- PRF set kept\n'); end
        cpf_report(sel);
    end
    prfUsed(scan) = prf;
    accIn = [0 0];  accOut = [0 0];  accNc = [0 0];
    tIn = zeros(1, nT);  tN = zeros(1, nT);  tTrue = cell(1, nT);  tFilt = zeros(1, nT);
    done = false;
    while ~done
        [blk, sim, truth] = csim_block(sim, prf);
        [y, st, info] = cpf_process(st, blk);
        if ~isempty(stNoComp)
            [y2, stNoComp] = cpf_process(stNoComp, blk);
        end
        % clutter region away from the targets
        cl = blk.rangeM > 1e3 & blk.rangeM < 19e3;
        tmask = false(size(blk.x));
        for j = 1:nT
            a = abs(mod(blk.azDeg - truth.tgt(j).azDeg + 180, 360) - 180) < 2 * S.beamwidthDeg;
            c = max(truth.tgt(j).cell - 4, 1):min(truth.tgt(j).cell + 4, S.nRange);
            tmask(a, c) = true;
            % target power after the filter: pulses within half a beam, minus the
            % power of the cells around it (noise + clutter residue)
            % (the target cell and its two neighbours: it moves during a block)
            on = abs(mod(blk.azDeg - truth.tgt(j).azDeg + 180, 360) - 180) < S.beamwidthDeg / 4;
            cj = truth.tgt(j).cell;
            if any(on) && cj > 10 && cj < S.nRange - 10
                ring = [cj-10:cj-4, cj+4:cj+10];
                tIn(j) = tIn(j) + sum(sum(abs(y(on, cj-1:cj+1)).^2));
                tN(j)  = tN(j) + 3 * sum(sum(abs(y(on, ring)).^2)) / numel(ring);
                tTrue{j} = truth.tgt(j);
                tFilt(j) = info.filter(find(on, 1), cj);
            end
        end
        m = ~tmask & repmat(cl, size(blk.x, 1), 1);
        accIn  = accIn  + [sum(abs(blk.x(m)).^2) nnz(m)];
        accOut = accOut + [sum(abs(y(m)).^2) nnz(m)];
        if ~isempty(stNoComp)
            accNc = accNc + [sum(abs(y2(m)).^2) nnz(m)];
        end
        done = truth.scanDone;
    end
    N0 = st.noise;
    resid(scan, :) = 10 * log10([accIn(1) / accIn(2), accOut(1) / accOut(2), accNc(1) / max(accNc(2), 1)] / N0);
    for j = 1:nT
        if tN(j) > 0
            tgtDb(scan, j) = 10 * log10(max(tIn(j) - tN(j), 1e-3 * tN(j)) / (tN(j) / 3));
            if tFilt(j) == 0
                g = 0;                                   % no clutter here: not filtered
            else
                bank = st.banks{find(abs(st.bankPrf - prf) < 1e-6, 1)};
                w = bank.W(:, tFilt(j));
                g = 10 * log10(abs(cpf_steering(2 * tTrue{j}.vRadMps / lam, prf, numel(w)) * w)^2);
            end
            fr = mod(tTrue{j}.rangeM / S.rangeCellM, 1);  % echo shared by two cells
            tgtPred(scan, j) = S.targets(j).snrDb + beamLossDb + 10 * log10((1 - fr)^2 + fr^2) + g;
        end
    end
    nc = '';
    if ~isempty(stNoComp), nc = sprintf('(%5.1f)', resid(scan, 3)); end
    fprintf('%3d  %5.0f  %4.1f m/s   %6.1f    %6.1f %-8s ', scan, prf, norm(sim.vel), ...
            resid(scan, 1), resid(scan, 2), nc);
    for j = 1:nT
        fprintf(' %5.1f(%5.1f)', tgtDb(scan, j), tgtPred(scan, j));
    end
    fprintf('\n');
end
fprintf('targets: ');
for j = 1:nT
    fprintf('%d: %.0f m/s  ', j, S.targets(j).vRadMps);
end
fprintf('\n');
m = st.map;
own = m.n >= cfg.map.minUpdates & m.C > st.noise;
fprintf('map: %d cells learned, CNR median %.1f dB, max %.1f dB, intrinsic sigma_v median %.2f m/s\n', ...
        nnz(own), median(10 * log10(m.C(own) / st.noise)), max(10 * log10(m.C(own) / st.noise)), ...
        median(sqrt(m.sig2(own & m.nw > 0))));

%% Plots
if doPlot
    figure('Name', 'PRF set');
    sel = st.sel;
    plot(sel.v, 10 * log10(sel.G.'));  hold on;
    plot(sel.v, 10 * log10(mean(sel.G, 1)), 'k', 'LineWidth', 2);
    grid on;  xlabel('radial velocity over ground [m/s]');  ylabel('filter gain [dB re noise]');
    legend([arrayfun(@(p) sprintf('%g Hz', p), sel.prfHz, 'UniformOutput', false), {'mean'}], ...
           'Location', 'southeast');
    ylim([-40 10]);
    title(sprintf('ripple of the mean %.1f dB, worst CA %.1f dB', sel.rippleDb, sel.caWorstDb));

    figure('Name', 'Clutter residue');
    plot(1:nScans, resid(:, 1), 'o-', 1:nScans, resid(:, 2), 's-', 1:nScans, resid(:, 3), 'x--');
    grid on;  xlabel('scan');  ylabel('dB re noise');
    legend('clutter in', 'after the filter', 'after the filter, no motion compensation');

    figure('Name', 'Clutter map');
    C = reshape(10 * log10(max(m.C, 1e-3 * st.noise) / st.noise), m.nA, m.nR);
    [ia, ir] = ndgrid((0.5:m.nA) * cfg.map.azCellDeg, (0.5:m.nR) * cfg.map.rangeCellM / 1e3);
    surf(ir .* sind(ia), ir .* cosd(ia), C, 'EdgeColor', 'none');  view(2);  axis equal;
    colorbar;  xlabel('east [km]');  ylabel('north [km]');  title('clutter map: CNR [dB]');
    xlim([-25 25]);  ylim([-25 25]);
end
