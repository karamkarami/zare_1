function scan = rsp_scan(P, S, varargin)
%RSP_SCAN Simulate and process antenna revolutions block by block, with a live PPI.
%
%   scan = rsp_scan(P, S, 'option', value, ...)
%
%   P : chain parameters (rsp_default_params)
%   S : scenario (see rsp_simulate)
%
%   Options
%     'degrees'  360     antenna rotation to process (720 = two scans)
%     'display'  true    live PPI while processing
%     'quiet'    false   no printing
%
%   Output
%     scan.det     every output detection (after the max stage):
%                  pulse, cell, bin, value, azDeg, rangeM, velocityMps
%     scan.plots   extracted plots (rsp_extract_plots), one per target echo
%     scan.truth   where each simulated target should appear
%     scan.report  truth against plots (matched plot and errors)
%     scan.maxValue, scan.maxBin   output stage as sparse pulses x range
%                  matrices (for rsp_ppi_replay)
%     scan.azDeg   azimuth of every pulse
%     scan.rangeM  range of every cell
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
R    = S.nRange;
g    = rsp_geometry(P, R);
P.output.keepDoppler = false;

nTotal = round(opt.degrees / g.degPerPulse);
hop    = P.fft.hop;
H      = numel(rsp_canceler_taps(P.canceler)) - 1 + P.fft.nPulses - 1 + ...
         (P.nci.nFrames - 1) * hop;                            % memory of the chain
H      = ceil(H / hop) * hop;
B      = max(hop, round(P.ppi.blockPulses / hop) * hop);
margin = ceil(1.5 * P.antenna.beamwidthDeg / g.degPerPulse);    % a plot is complete after this

det  = struct('pulse', zeros(0, 1), 'cell', zeros(0, 1), 'bin', zeros(0, 1), 'value', zeros(0, 1));
hist = [];  histIdx = zeros(0, 1);
sim  = [];  cache = [];
lastEmit = 0;
if opt.display
    ppi = rsp_ppi_init(P, g.rangeM, 'PPI');
end
if ~opt.quiet
    fprintf('Scan: %d pulses (%.0f deg), %d per block, %.3f deg per pulse\n', ...
            nTotal, opt.degrees, B, g.degPerPulse);
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
    v    = [hist; vNew];
    vIdx = [histIdx; pulses];
    o    = rsp_chain(v, P, cache);
    cache = o.cache;

    % rows of this block's new pulses
    pr  = vIdx(o.idx.max);
    new = pr >= max(1, pulses(1));
    pr  = pr(new);
    val = o.max.value(new, :);
    bin = o.max.bin(new, :);
    [ii, cc] = find(val > 0);
    lin = sub2ind(size(val), ii, cc);
    det.pulse = [det.pulse; pr(ii)];
    det.cell  = [det.cell;  cc];
    det.bin   = [det.bin;   double(bin(lin))];
    det.value = [det.value; double(val(lin))];

    % plots completed in this block (for the display)
    if opt.display
        final = next > nTotal;
        upto  = pr(end) - margin * ~final;
        win   = det.pulse > lastEmit - 4*margin;
        pl    = rsp_extract_plots(subset(det, win), P, R);
        emit  = pl.pulse > lastEmit & pl.pulse <= upto;
        pl    = subset(pl, emit);
        lastEmit = max(lastEmit, upto);

        vr = o.idx.max(new);                                    % rows of v
        switch lower(P.ppi.source)
            case 'max',      rows = 10*log10(val);
            case 'video',    rows = 20*log10(abs(double(v(vr, :))));
            case 'decoder',  rows = 20*log10(abs(double(o.decoder(vr, :))));
            case 'canceler'
                [~, cr] = ismember(vr, o.idx.canceler);
                rows = 20*log10(abs(double(o.canceler(cr, :))));
            otherwise, error('rsp_scan:source', 'Unknown P.ppi.source "%s".', P.ppi.source);
        end
        gp  = rsp_geometry(P, R, pr);
        ppi = rsp_ppi_update(ppi, gp.azDeg, rows, pl, ...
              sprintf('Az %5.1f deg   t = %5.1f s   %s   bins %s', gp.azDeg(end), ...
                      gp.timeS(end), upper(P.ppi.source), binText(P.output.bins)));
    end

    hist    = v(end-H+1:end, :);
    histIdx = vIdx(end-H+1:end);
    if ~opt.quiet
        fprintf('  pulses %6d .. %6d  az %6.1f deg  %6d detections  (%.1f s)\n', ...
                max(pulses(1), 1), pulses(end), mod(g.degPerPulse*(pulses(end)-1) + ...
                P.radar.azStartDeg, 360), numel(det.pulse), toc(t0));
    end
end

% --- results -------------------------------------------------------------------------------
gAll        = rsp_geometry(P, R, (1:nTotal)');
scan.azDeg  = gAll.azDeg;
scan.rangeM = g.rangeM;
det.azDeg   = gAll.azDeg(det.pulse);
det.rangeM  = g.rangeM(det.cell);
det.velocityMps = g.binVelMps(det.bin)';
scan.det    = det;
scan.maxValue = sparse(det.pulse, det.cell, det.value, nTotal, R);
scan.maxBin   = sparse(det.pulse, det.cell, det.bin,   nTotal, R);
scan.plots  = rsp_extract_plots(det, P, R);
scan.truth  = rsp_truth(P, S);
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
    printReport(scan.report);
end
end

% ==========================================================================================
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
% Nearest plot to every target (within 2 beamwidths and 500 m).
rep = struct('target', {}, 'found', {}, 'azDeg', {}, 'rangeM', {}, 'velocityMps', {}, ...
             'plotAzDeg', {}, 'plotRangeM', {}, 'plotVelocityMps', {}, ...
             'dAzDeg', {}, 'dRangeM', {}, 'dVelMps', {});
for k = 1:numel(truth)
    t  = truth(k);
    da = mod(plots.azDeg - t.azDeg + 180, 360) - 180;
    dr = plots.rangeM - t.rangeM;
    ok = abs(da) < 2*P.antenna.beamwidthDeg & abs(dr) < 500;
    rep(k).target      = k;
    rep(k).found       = any(ok);
    rep(k).azDeg       = t.azDeg;
    rep(k).rangeM      = t.rangeM;
    rep(k).velocityMps = t.velocityMps;
    if any(ok)
        cand = find(ok);
        [~, j] = min(abs(da(cand)) / P.antenna.beamwidthDeg + abs(dr(cand)) / 500);
        j = cand(j);
        rep(k).plotAzDeg       = plots.azDeg(j);
        rep(k).plotRangeM      = plots.rangeM(j);
        rep(k).plotVelocityMps = plots.velocityMps(j);
        rep(k).dAzDeg  = da(j);
        rep(k).dRangeM = dr(j);
        dv = plots.velocityMps(j) - t.foldedMps;
        rep(k).dVelMps = mod(dv + g.vUnambMps/2, g.vUnambMps) - g.vUnambMps/2;
    end
end
end

function printReport(rep)
fprintf('\n  target |  az [deg]  range [m]  vel [m/s] |   plot az   plot range  plot vel | dAz [deg]  dR [m]  dV [m/s]\n');
for k = 1:numel(rep)
    r = rep(k);
    if r.found
        fprintf('  %6d | %8.2f %10.0f %9.1f | %9.2f %11.0f %9.1f | %8.2f %7.0f %8.1f\n', ...
            r.target, r.azDeg, r.rangeM, r.velocityMps, r.plotAzDeg, r.plotRangeM, ...
            r.plotVelocityMps, r.dAzDeg, r.dRangeM, r.dVelMps);
    else
        fprintf('  %6d | %8.2f %10.0f %9.1f |   not detected\n', r.target, r.azDeg, ...
            r.rangeM, r.velocityMps);
    end
end
fprintf('  (plot velocity is folded into +-lambda*PRF/4; dV compares with the folded truth)\n\n');
end
