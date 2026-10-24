function fa = rsp_false_alarms(scan, show)
%RSP_FALSE_ALARMS False alarms of a scan, per clutter region, against the design Pfa.
%
%   fa = rsp_false_alarms(scan)
%   fa = rsp_false_alarms(scan, true)     also print the table
%
%   scan : result of rsp_scan
%
%   Detections (output after the max over the Doppler bins) that are not
%   near a target crossing are false detections. A gate of half the
%   budgeted target width + one beamwidth in azimuth and 1 km in range
%   (apparent range) around each crossing is left out. Plots that match no
%   target are false plots. Both are counted per clutter region
%   (rsp_clutter_region, regions widened by one beamwidth and 1 km).
%
%   fa.names         region names: 'clear', the layers, 'shadow'
%   fa.cells         output cells tested in each region (rows x range cells)
%   fa.falseDet      false detections in each region
%   fa.pfa           measured false alarm probability per output cell
%   fa.pfaDesign     design value per output cell: 1 - (1 - P.cfar.pfa)^nBins,
%                    nBins = bins that enter the max (NaN with a fixed factor)
%   fa.falsePlots    false plots in each region,  fa.plotsPerScan  per revolution
%   fa.plots         the false plots (azDeg, rangeM, velocityMps, marginDb, region)
%   fa.isFalse       logical, one per scan.det entry
%   fa.found, fa.echoes, fa.pd   target echoes found / expected, ratio
%   fa.nScans        revolutions processed
%   (the totals are in the last element of cells / falseDet / pfa /
%    falsePlots, named 'total')

if nargin < 2
    show = false;
end
P   = scan.P;
R   = numel(scan.rangeM);
g   = rsp_geometry(P, R);
bw  = P.antenna.beamwidthDeg;
L   = scan.clutter;
det = scan.det;
dly = g.outputDelayPulses;
nb  = numel(P.output.bins);
if nb == 0
    nb = max(P.fft.nfft, P.fft.nPulses);
end

% --- detections near a target crossing are not false -------------------------------
azC = mod(det.azDeg - g.degPerPulse * dly, 360);               % centre of the data
pC  = det.pulse - dly;
isF = true(size(det.pulse));
for k = 1:numel(scan.truth)
    t  = scan.truth(k);
    w  = bw;
    if ~isempty(scan.budget) && t.target <= numel(scan.budget)
        w = max(scan.budget(t.target).widthDeg / 2, bw) + bw;
    end
    da = abs(mod(azC - t.azDeg + 180, 360) - 180);
    isF = isF & ~(da < w & abs(det.rangeM - t.apparentRangeM) < 1000 & ...
                  abs(pC - t.pulse) < 2 * w / g.degPerPulse);
end

% --- regions ------------------------------------------------------------------------
[idDet, names] = rsp_clutter_region(L, azC(isF), det.rangeM(isF), bw, 1000);
nR = numel(names);
falseDet = count(idDet, nR);

% tested cells: every output row x every CFAR cell, mapped on a 0.5 deg grid
cells = scan.validCells(1):scan.validCells(2);
azAll = mod(scan.azDeg - g.degPerPulse * dly, 360);
nAz   = 720;
hc    = accumarray(min(floor(azAll / (360/nAz)), nAz - 1) + 1, 1, [nAz 1]);
azG   = ((0:nAz-1)' + 0.5) * 360/nAz;
used  = find(hc > 0);
nCell = zeros(nR, 1);
for i = 1:numel(used)
    id = rsp_clutter_region(L, repmat(azG(used(i)), 1, numel(cells)), ...
                            scan.rangeM(cells).', bw, 1000);
    nCell = nCell + hc(used(i)) * count(id, nR);
end

% --- false plots ----------------------------------------------------------------------
pl   = scan.plots;
used = false(size(pl.azDeg));
rep  = scan.report;
for k = 1:numel(rep)
    if rep(k).found
        used(rep(k).plot) = true;
    end
end
fp = find(~used);
plAz = mod(pl.azDeg(fp), 360);
idPl = rsp_clutter_region(L, plAz, pl.rangeM(fp), bw, 1000);
falsePlots = count(idPl, nR);

fa.names      = [names, {'total'}];
fa.cells      = [nCell; sum(nCell)];
fa.falseDet   = [falseDet; sum(falseDet)];
fa.pfa        = fa.falseDet ./ max(fa.cells, 1);
if strcmpi(P.cfar.thresholdMode, 'pfa')
    fa.pfaDesign = 1 - (1 - P.cfar.pfa)^nb;
else
    fa.pfaDesign = NaN;
end
fa.nScans       = numel(scan.azDeg) * g.degPerPulse / 360;
fa.falsePlots   = [falsePlots; sum(falsePlots)];
fa.plotsPerScan = fa.falsePlots / max(fa.nScans, eps);
fa.plots = struct('azDeg', plAz, 'rangeM', pl.rangeM(fp), 'velocityMps', pl.velocityMps(fp), ...
                  'marginDb', pl.marginDb(fp), 'region', idPl(:));
fa.isFalse = isF;
fa.found   = nnz([rep.found]);
fa.echoes  = numel(rep);
fa.pd      = fa.found / max(fa.echoes, 1);

if show
    fprintf(['False alarms (output after the max over %d bins; design Pfa %.1e per bin ' ...
             '-> %.1e per output cell)\n'], nb, P.cfar.pfa, fa.pfaDesign);
    fprintf('  %-12s %14s %10s %13s %12s %10s\n', 'region', 'tested cells', 'false det', ...
            'Pfa measured', 'false plots', 'per scan');
    for i = 1:numel(fa.names)
        if fa.cells(i) == 0 && fa.falseDet(i) == 0 && i < numel(fa.names)
            continue
        end
        fprintf('  %-12s %14d %10d %13.2e %12d %10.2f\n', fa.names{i}, fa.cells(i), ...
                fa.falseDet(i), fa.pfa(i), fa.falsePlots(i), fa.plotsPerScan(i));
    end
    fprintf('  Target echoes detected: %d of %d (Pd %.2f) in %.2f revolution(s)\n\n', ...
            fa.found, fa.echoes, fa.pd, fa.nScans);
end
end

% ======================================================================================
function n = count(id, nR)
% occurrences of the region ids 0 .. nR-1 (column)
n = sum(bsxfun(@eq, double(id(:)), 0:nR-1), 1)';
end
