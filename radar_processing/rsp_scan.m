function scan = rsp_scan(P, S, varargin)
%RSP_SCAN Rotating radar: simulate and process revolutions, with a live PPI.
%
%   scan = rsp_scan(P, S, 'option', value, ...)
%
%   P : chain parameters (radar_params / rsp_default_params)
%   S : scenario (see rsp_simulate)
%
%   Options
%     'degrees'  360     antenna rotation to process (720 = two revolutions)
%     'display'  true    live PPI: the sweep turns and paints as the pulses
%                        are processed (P.ppi.stepDeg steps, at P.ppi.speed
%                        times the real antenna speed when processing keeps up)
%     'quiet'    false   no printing
%
%   Output
%     scan.det     every output detection (after the max stage): pulse, cell,
%                  bin, binFrac, value, marginDb, azDeg, rangeM, velocityMps
%     scan.plots   extracted plots (rsp_extract_plots), one per target echo
%     scan.truth   where each target should appear, every revolution
%     scan.report  truth against plots: found, errors, CFAR margin
%     scan.budget  expected margins (rsp_budget)
%     scan.maxValue, scan.maxMargin   output stage as sparse pulses x range
%                  matrices (for rsp_ppi_replay)
%     scan.azDeg   azimuth of every pulse;  scan.rangeM  range of every cell
%     scan.ppi     PPI state (if displayed)
%
%   The video is generated and processed in blocks of P.ppi.blockPulses
%   pulses. Each block is preceded by the last pulses of the previous one
%   (canceler + FFT + integration memory), so the output is identical to
%   processing the whole scan at once, with bounded memory.

opt = struct('degrees', 360, 'display', true, 'quiet', false);
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
R = S.nRange;
g = rsp_geometry(P, R);
P.output.keepDoppler = false;

nTotal = round(opt.degrees / g.degPerPulse);
nScans = ceil(nTotal / (360 / g.degPerPulse) - 1e-9);
hop    = P.fft.hop;
H      = numel(rsp_canceler_taps(P.canceler)) - 1 + P.fft.nPulses - 1 + ...
         (P.nci.nFrames - 1) * hop;                            % memory of the chain
H      = ceil(H / hop) * hop;
B      = max(hop, round(P.ppi.blockPulses / hop) * hop);
lag    = ceil(1.5 * P.antenna.beamwidthDeg / g.degPerPulse);   % a plot is complete after this
stepP  = max(1, round(P.ppi.stepDeg / g.degPerPulse));         % pulses per screen update

empty = zeros(0, 1);
det  = struct('pulse', empty, 'cell', empty, 'bin', empty, 'binFrac', empty, ...
              'value', empty, 'marginDb', empty);
hist = [];  histIdx = empty;
sim  = [];  cache = [];
lastEmit = 0;

scan.budget = rsp_budget(P, S, opt.quiet);
if opt.display
    ppi = rsp_ppi_init(P, g.rangeM, 'PPI');
end
if ~opt.quiet
    fprintf('Scan: %d pulses (%.0f deg, %d revolution(s)), %d per block, %.3f deg per pulse\n', ...
            nTotal, opt.degrees, nScans, B, g.degPerPulse);
end

t0   = tic;
next = 1 - H;                                                   % pre-roll fills the memory
while next <= nTotal
    if isempty(hist)
        pulses = (next:min(B, nTotal))';
    else
        pulses = (next:min(next + B - 1, nTotal))';
    end
    next = pulses(end) + 1;

    [vNew, sim] = rsp_simulate(P, S, pulses, sim);
    v     = [hist; vNew];
    vIdx  = [histIdx; pulses];
    o     = rsp_chain(v, P, cache);
    cache = o.cache;

    % output rows of this block's new pulses
    pr  = vIdx(o.idx.max);
    new = pr >= max(1, pulses(1));
    pr  = pr(new);
    val = o.max.value(new, :);
    [ii, cc] = find(val > 0);
    lin = sub2ind(size(val), ii, cc);
    bin = o.max.bin(new, :);       bf = o.max.binFrac(new, :);   mg = o.max.marginDb(new, :);
    det.pulse    = [det.pulse;    pr(ii)];
    det.cell     = [det.cell;     cc];
    det.bin      = [det.bin;      double(bin(lin))];
    det.binFrac  = [det.binFrac;  double(bf(lin))];
    det.value    = [det.value;    double(val(lin))];
    det.marginDb = [det.marginDb; double(mg(lin))];

    if opt.display
        % plots completed by the end of this block
        final = next > nTotal;
        upto  = pr(end) - lag * ~final;
        win   = det.pulse > lastEmit - 4*lag;
        pl    = rsp_extract_plots(subset(det, win), P, R);
        pl    = subset(pl, pl.pulse > lastEmit & pl.pulse <= upto);
        lastEmit = max(lastEmit, upto);

        vr = o.idx.max(new);                                    % rows of v
        switch lower(P.ppi.source)
            case 'max'
                if strcmpi(P.ppi.maxScale, 'margin')
                    rows = double(mg);                          % dB over the threshold
                    rows(~isfinite(rows)) = -Inf;
                else
                    rows = 10*log10(double(val));
                end
            case 'video',    rows = 20*log10(abs(double(v(vr, :))));
            case 'decoder',  rows = 20*log10(abs(double(o.decoder(vr, :))));
            case 'canceler'
                [~, cr] = ismember(vr, o.idx.canceler);
                rows = 20*log10(abs(double(o.canceler(cr, :))));
            otherwise, error('rsp_scan:source', 'Unknown P.ppi.source "%s".', P.ppi.source);
        end
        gp = rsp_geometry(P, R, pr);

        % paint the block in small steps: the sweep turns smoothly, paced to
        % the antenna speed when the processing is fast enough
        for s0 = 1:stepP:numel(pr)
            k  = (s0:min(s0 + stepP - 1, numel(pr)))';
            sp = pl.pulse + lag <= pr(k(end)) | (k(end) == numel(pr));
            ppi = rsp_ppi_update(ppi, gp.azDeg(k), rows(k, :), subset(pl, sp), ...
                  sprintf('Az %5.1f deg   t = %5.1f s   scan %d   %s   bins %s', ...
                          gp.azDeg(k(end)), gp.timeS(k(end)), ...
                          floor((pr(k(end)) - 1) * g.degPerPulse / 360) + 1, ...
                          upper(P.ppi.source), binText(P.output.bins)));
            pl  = subset(pl, ~sp);
            wait = gp.timeS(k(end)) / P.ppi.speed - toc(t0);
            if isfinite(P.ppi.speed) && wait > 0
                pause(wait);
            end
        end
    end

    hist    = v(end-H+1:end, :);
    histIdx = vIdx(end-H+1:end);
    if ~opt.quiet
        fprintf('  pulses %6d .. %6d  az %6.1f deg  %6d detections  (%.1f s)\n', ...
                max(pulses(1), 1), pulses(end), mod(g.degPerPulse*(pulses(end)-1) + ...
                P.radar.azStartDeg, 360), numel(det.pulse), toc(t0));
    end
end

% --- results ----------------------------------------------------------------------------
gAll        = rsp_geometry(P, R, (1:nTotal)');
scan.azDeg  = gAll.azDeg;
scan.rangeM = g.rangeM;
det.azDeg   = gAll.azDeg(det.pulse);
det.rangeM  = g.rangeM(det.cell);
det.velocityMps = g.binVelMps(det.bin)';
scan.det    = det;
scan.maxValue  = sparse(det.pulse, det.cell, det.value, nTotal, R);
scan.maxMargin = sparse(det.pulse, det.cell, max(det.marginDb, 1e-3), nTotal, R);
scan.plots  = rsp_extract_plots(det, P, R);
tr          = rsp_truth(P, S, nScans);
scan.truth  = tr([tr.pulse] <= nTotal);
scan.report = matchTruth(scan.truth, scan.plots, P, g);
scan.P      = P;
scan.S      = S;
scan.time   = toc(t0);
if opt.display
    scan.ppi = ppi;
end
if ~opt.quiet
    fprintf('Done in %.1f s: %d detections, %d plots\n', scan.time, numel(det.pulse), ...
            numel(scan.plots.azDeg));
    printReport(scan.report, scan.plots, g);
end
end

% ========================================================================================
function s = subset(s, k)
% Rows k of every field of a struct of column vectors.
f = fieldnames(s);
for i = 1:numel(f)
    s.(f{i}) = s.(f{i})(k);
end
end

function t = binText(b)
if isempty(b)
    t = 'all';
else
    t = mat2str(b);
end
end

function rep = matchTruth(truth, plots, P, g)
% Nearest plot to every target crossing (within 2 beamwidths and 500 m).
rep = struct('target', {}, 'scan', {}, 'found', {}, 'azDeg', {}, 'rangeM', {}, ...
             'velocityMps', {}, 'plot', {}, 'plotAzDeg', {}, 'plotRangeM', {}, ...
             'plotVelocityMps', {}, 'dAzDeg', {}, 'dRangeM', {}, 'dVelMps', {}, ...
             'hits', {}, 'marginDb', {});
used = false(size(plots.azDeg));
for k = 1:numel(truth)
    t  = truth(k);
    da = mod(plots.azDeg - t.azDeg + 180, 360) - 180;
    dr = plots.rangeM - t.rangeM;
    dp = abs(plots.pulse - t.pulse);                             % same revolution
    ok = abs(da) < 2*P.antenna.beamwidthDeg & abs(dr) < 500 & dp < 90 / g.degPerPulse & ~used;
    rep(k).target      = t.target;
    rep(k).scan        = t.scan;
    rep(k).found       = any(ok);
    rep(k).azDeg       = t.azDeg;
    rep(k).rangeM      = t.rangeM;
    rep(k).velocityMps = t.velocityMps;
    if any(ok)
        cand = find(ok);
        [~, j] = min(abs(da(cand)) / P.antenna.beamwidthDeg + abs(dr(cand)) / 500);
        j = cand(j);
        used(j) = true;
        rep(k).plot            = j;
        rep(k).plotAzDeg       = plots.azDeg(j);
        rep(k).plotRangeM      = plots.rangeM(j);
        rep(k).plotVelocityMps = plots.velocityMps(j);
        rep(k).dAzDeg   = da(j);
        rep(k).dRangeM  = dr(j);
        dv = plots.velocityMps(j) - t.foldedMps;
        rep(k).dVelMps  = mod(dv + g.vUnambMps/2, g.vUnambMps) - g.vUnambMps/2;
        rep(k).hits     = plots.hits(j);
        rep(k).marginDb = plots.marginDb(j);
    end
end
end

function printReport(rep, plots, g)
fprintf(['\n  target scan |  az[deg]  range[m]  vel[m/s] |  plot az  plot range  plot vel | ' ...
         ' dAz[deg]  dR[m]  dV[m/s] | hits  CFAR margin[dB]\n']);
for k = 1:numel(rep)
    r = rep(k);
    if r.found
        fprintf('  %6d %4d | %8.2f %9.0f %9.1f | %8.2f %11.0f %9.1f | %8.2f %6.0f %8.1f | %4d %10.1f\n', ...
            r.target, r.scan, r.azDeg, r.rangeM, r.velocityMps, r.plotAzDeg, r.plotRangeM, ...
            r.plotVelocityMps, r.dAzDeg, r.dRangeM, r.dVelMps, r.hits, r.marginDb);
    else
        fprintf('  %6d %4d | %8.2f %9.0f %9.1f |   NOT DETECTED\n', r.target, r.scan, ...
            r.azDeg, r.rangeM, r.velocityMps);
    end
end
nFalse = numel(plots.azDeg) - nnz([rep.found]);
fprintf(['  %d of %d target echoes detected, %d other plot(s). Velocities are measured ' ...
         'within +-%.0f m/s.\n\n'], nnz([rep.found]), numel(rep), nFalse, g.vUnambMps/2);
end
