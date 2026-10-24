function files = rsp_save(res, folder, varargin)
%RSP_SAVE Save the results of a scan or a schedule: .mat, CSV tables, text and figures.
%
%   files = rsp_save(res, folder, 'option', value, ...)
%
%   res    : result of rsp_scan (one scan) or rsp_schedule (run.scans)
%   folder : output folder; '' = results/run_yyyymmdd_HHMMSS next to this file
%
%   Options
%     'figures'     true   save the PPI (if displayed) and the summary figures
%     'detections'  true   write every detection to detections.csv (large)
%
%   Files
%     results.mat        everything (without graphics handles): variable 'result'
%     plots.csv          entry, azDeg, rangeM, velocityMps, marginDb, powerDb, hits,
%                        widthDeg, isFalse
%     detections.csv     entry, pulse, cell, azDeg, rangeM, bin, velocityMps,
%                        value, marginDb, isFalse
%     report.csv         every target echo against its plot (errors, margin)
%     budget.csv         link budget of every target
%     falsealarms.csv    per entry and clutter region: tested cells, false
%                        detections, measured and design Pfa, false plots
%     timeline.csv       schedule timeline (rsp_schedule only)
%     summary.txt        the printed tables
%     ppi.png, summary_<k>.png, scenario.png
%
%   Read back:  s = load(fullfile(folder, 'results.mat'));  s.result
%               t = dlmread(fullfile(folder, 'plots.csv'), ',', 1, 0);

opt = struct('figures', true, 'detections', true);
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
if nargin < 2 || isempty(folder)
    folder = fullfile(fileparts(mfilename('fullpath')), 'results', ...
                      ['run_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
if ~exist(folder, 'dir')
    mkdir(folder);
end
files = {};

isRun = isfield(res, 'scans');
if isRun
    scans = res.scans;
else
    scans = res;
end
ppi = [];
if isfield(res, 'ppi')
    ppi = res.ppi;
end

% --- .mat ------------------------------------------------------------------------------------
result = res;
if isfield(result, 'ppi')
    result = rmfield(result, 'ppi');
end
if isRun
    result.sched = rmHandles(result.sched);
end
f = fullfile(folder, 'results.mat');
save(f, 'result', '-v7');
files{end+1} = f;

% --- CSV tables -----------------------------------------------------------------------------
pl = [];  dt = [];  rp = [];  bu = [];
faLines = {};
for k = 1:numel(scans)
    s  = scans(k);
    fa = s.falseAlarms;
    p  = s.plots;
    isF = true(size(p.azDeg));
    for j = 1:numel(s.report)
        if s.report(j).found
            isF(s.report(j).plot) = false;
        end
    end
    pl = [pl; repmat(k, numel(p.azDeg), 1), p.azDeg, p.rangeM, p.velocityMps, p.marginDb, ...
          p.powerDb, p.hits, p.widthDeg, isF];                               %#ok<AGROW>
    if opt.detections
        d  = s.det;
        dt = [dt; repmat(k, numel(d.pulse), 1), d.pulse, d.cell, d.azDeg, d.rangeM, d.bin, ...
              d.velocityMps(:), d.value, d.marginDb, fa.isFalse];          %#ok<AGROW>
    end
    for j = 1:numel(s.report)
        r = s.report(j);
        if r.found
            v = [r.plotAzDeg r.plotRangeM r.plotVelocityMps r.dAzDeg r.dRangeM r.dVelMps ...
                 r.hits r.marginDb];
        else
            v = nan(1, 8);
        end
        rp = [rp; k r.target r.scan r.found r.azDeg r.rangeM r.velocityMps v r.shadowDb];  %#ok<AGROW>
    end
    for j = 1:numel(s.budget)
        b = s.budget(j);
        bu = [bu; k b.target b.rangeM b.velocityMps b.snrInDb b.pulse b.decoderDb b.mtiDb ...
              b.fftDb b.snrOutDb b.bin b.thresholdDb b.marginDb b.widthDeg b.foldedMps];  %#ok<AGROW>
    end
    for j = 1:numel(fa.names)
        faLines{end+1} = sprintf('%d,%s,%d,%d,%.6g,%.6g,%d,%.6g', k, fa.names{j}, ...
                                 fa.cells(j), fa.falseDet(j), fa.pfa(j), fa.pfaDesign, ...
                                 fa.falsePlots(j), fa.plotsPerScan(j));        %#ok<AGROW>
    end
end
files{end+1} = writeCsv(fullfile(folder, 'plots.csv'), ['entry,azDeg,rangeM,velocityMps,' ...
               'marginDb,powerDb,hits,widthDeg,isFalse'], pl);
if opt.detections
    files{end+1} = writeCsv(fullfile(folder, 'detections.csv'), ['entry,pulse,cell,azDeg,' ...
                   'rangeM,bin,velocityMps,value,marginDb,isFalse'], dt);
end
files{end+1} = writeCsv(fullfile(folder, 'report.csv'), ['entry,target,scan,found,azDeg,' ...
               'rangeM,velocityMps,plotAzDeg,plotRangeM,plotVelocityMps,dAzDeg,dRangeM,' ...
               'dVelMps,hits,marginDb,shadowDb'], rp);
files{end+1} = writeCsv(fullfile(folder, 'budget.csv'), ['entry,target,rangeM,velocityMps,' ...
               'snrInDb,pulse,decoderDb,mtiDb,fftDb,snrOutDb,bin,thresholdDb,marginDb,' ...
               'widthDeg,foldedMps'], bu);
files{end+1} = writeLines(fullfile(folder, 'falsealarms.csv'), ...
               [{'entry,region,cells,falseDet,pfa,pfaDesign,falsePlots,plotsPerScan'}, faLines]);
if isRun
    tl = res.timeline;
    lines = {['entry,name,startS,endS,durationS,startAzDeg,endAzDeg,turns,prfHz,fcMHz,' ...
              'lambdaM,rangeUnambKm,vUnambMps,canceler,fftPulses,nciFrames,pulses,found,' ...
              'echoes,falsePlots,pfa,pfaDesign']};
    for k = 1:numel(tl.entry)
        lines{end+1} = sprintf(['%d,%s,%.4f,%.4f,%.4f,%.3f,%.3f,%.4f,%g,%g,%.5f,%.3f,%.3f,' ...
                                '%d,%d,%d,%d,%d,%d,%d,%.6g,%.6g'], tl.entry(k), ...
                               strrep(tl.name{k}, ',', ';'), tl.startS(k), tl.endS(k), ...
                               tl.durationS(k), tl.startAzDeg(k), tl.endAzDeg(k), tl.turns(k), ...
                               tl.prfHz(k), tl.fcMHz(k), tl.lambdaM(k), tl.rangeUnambKm(k), ...
                               tl.vUnambMps(k), tl.canceler(k), tl.fftPulses(k), ...
                               tl.nciFrames(k), tl.pulses(k), tl.found(k), tl.echoes(k), ...
                               tl.falsePlots(k), tl.pfa(k), tl.pfaDesign(k));   %#ok<AGROW>
    end
    files{end+1} = writeLines(fullfile(folder, 'timeline.csv'), lines);
end

% --- summary text ---------------------------------------------------------------------------
txt = '';
if isRun
    txt = evalc('rsp_schedule_print(res.timeline)');
end
for k = 1:numel(scans)
    nm = sprintf('entry %d', k);
    if isfield(scans(k), 'name') && ~isempty(scans(k).name)
        nm = scans(k).name;
    end
    txt = [txt sprintf('\n--- %d: %s ---\n', k, nm) ...
           evalc('rsp_false_alarms(scans(k), true);')];                 %#ok<AGROW>
end
files{end+1} = writeLines(fullfile(folder, 'summary.txt'), {txt});

% --- figures ------------------------------------------------------------------------------
if opt.figures
    if ~isempty(ppi) && isfield(ppi, 'fig') && ishandle(ppi.fig)
        f = fullfile(folder, 'ppi.png');
        print(ppi.fig, f, '-dpng', '-r100');
        files{end+1} = f;
    end
    for k = 1:numel(scans)
        fig = rsp_plot_summary(scans(k), 'visible', 'off');
        f = fullfile(folder, sprintf('summary_%d.png', k));
        print(fig, f, '-dpng', '-r100');
        close(fig);
        files{end+1} = f;                                                  %#ok<AGROW>
    end
    fig = rsp_plot_scenario(scans(1).P, scans(1).S, 'visible', 'off');
    f = fullfile(folder, 'scenario.png');
    print(fig, f, '-dpng', '-r100');
    close(fig);
    files{end+1} = f;
end
fprintf('Saved %d files to %s\n', numel(files), folder);
end

% ======================================================================================
function f = writeCsv(f, header, M)
fid = fopen(f, 'w');
fprintf(fid, '%s\n', header);
if ~isempty(M)
    fmt = [repmat('%.10g,', 1, size(M, 2) - 1) '%.10g\n'];
    fprintf(fid, fmt, double(M).');
end
fclose(fid);
end

function f = writeLines(f, lines)
fid = fopen(f, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end

function s = rmHandles(s)
% function handles (sched.apply) are saved as text
if isfield(s, 'apply')
    for k = 1:numel(s)
        if isa(s(k).apply, 'function_handle')
            s(k).apply = func2str(s(k).apply);
        end
    end
end
end
