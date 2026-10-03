function rsp_plot(out, video, varargin)
%RSP_PLOT Show the output of every block of the chain.
%
%   rsp_plot(out, video, 'option', value, ...)
%
%   out   : result of rsp_chain
%   video : chain input ([] = not shown)
%
%   Options
%     'frame'    []      row of the integral / CFAR lanes for the range cut
%                        ([] = row with the most detections)
%     'bin'      []      Doppler bin of the integral panel and the range cut
%                        ([] = bin of the strongest detection)
%     'bins'     []      bins of the output panel (max after CFAR);
%                        [] = as processed (P.output.bins)
%     'range'    []      range cells shown, [first last] ([] = all)
%     'colormap' 'jet'

opt = struct('frame', [], 'bin', [], 'bins', [], 'range', [], 'colormap', 'jet');
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
R = size(out.decoder, 2);
if isempty(opt.range)
    rc = 1:R;
else
    rc = opt.range(1):min(opt.range(end), R);
end

det = out.cfar.det;
if isempty(opt.frame)
    [~, opt.frame] = max(sum(sum(det, 3), 2));
end
if isempty(opt.bin)
    d = squeeze(det(opt.frame, :, :));
    s = squeeze(out.integral(opt.frame, :, :));
    s(~d) = -Inf;
    [~, k] = max(s(:));
    [~, opt.bin] = ind2sub(size(s), k);
end

toDb = @(x) 20*log10(max(abs(double(x)), realmin));
pwDb = powerDb(out.P.nci.law);

figure('Name', 'Radar processing chain', 'Color', 'w', 'Position', [60 60 1500 800]);
panel = 0;
if ~isempty(video)
    panel = panel + 1;
    image2(panel, toDb(video(:, rc)), rc, 'video');
end
image2(panel + 1, toDb(out.decoder(:, rc)), rc, sprintf('matched filter, switch %d', ...
       out.mf.switchCell));
image2(panel + 2, toDb(out.canceler(:, rc)), rc, sprintf('%d-pulse canceler', ...
       numel(rsp_canceler_taps(out.P.canceler))));
image2(panel + 3, pwDb(out.integral(:, rc, opt.bin)), rc, ...
       sprintf('integral n=%d, bin %d', out.P.nci.nFrames, opt.bin));

% output stage: max over the selected bins after the CFAR
if isempty(opt.bins)
    val = out.max.value;
    bt  = 'all';
    if ~isempty(out.P.output.bins), bt = mat2str(out.P.output.bins); end
else
    val = rsp_cfar_max(out.cfar.map, opt.bins);
    bt  = mat2str(opt.bins);
end
subplot(2, 3, panel + 4);
[fr, cc] = find(val(:, rc) > 0);
lvl = pwDb(val(sub2ind(size(val), fr, cc + rc(1) - 1)));
scatter(cc + rc(1) - 1, fr, 8, lvl, 'filled');
axis([rc(1) rc(end) 0.5 size(val, 1) + 0.5]);
box on; cb = colorbar; ylabel(cb, 'dB');
xlabel('range cell'); ylabel('row');
title(sprintf('%s-CFAR + max, bins %s: %d', upper(out.P.cfar.type), bt, nnz(val(:, rc) > 0)));

% range cut: integrated signal against the CFAR threshold
subplot(2, 3, panel + 5);
x = squeeze(out.integral(opt.frame, rc, opt.bin));
t = squeeze(out.cfar.threshold(opt.frame, rc, opt.bin));
h = squeeze(det(opt.frame, rc, opt.bin));
plot(rc, pwDb(x), 'Color', [0.2 0.4 0.8]); hold on
plot(rc, pwDb(t), 'r', 'LineWidth', 1);
plot(rc(h), pwDb(x(h)), 'ko', 'MarkerFaceColor', 'y');
hold off; grid on; xlim([rc(1) rc(end)]);
v = pwDb(x);  v = sort(v(isfinite(v)));
if ~isempty(v)
    tt = pwDb(t(isfinite(t)));
    ylim([v(max(1, round(0.01*numel(v)))) - 5, max([v(end); tt(:)]) + 5]);
end
xlabel('range cell'); ylabel('dB');
legend('integral', 'threshold', 'detection', 'Location', 'best');
title(sprintf('frame %d (pulse %d), bin %d', opt.frame, ...
      out.idx.integral(opt.frame), opt.bin));
colormap(opt.colormap);
end

function image2(k, img, rc, ttl)
subplot(2, 3, k);
imagesc(rc, 1:size(img, 1), img, rsp_db_limits(img));
axis xy; colorbar;
xlabel('range cell'); ylabel('row');
title([ttl ' [dB]']);
end

function f = powerDb(law)
% dB conversion matching the detector law of the integrator.
switch lower(law)
    case 'square', f = @(x) 10*log10(nanZero(x));
    case 'linear', f = @(x) 20*log10(nanZero(x));
    otherwise,     f = @(x) double(x);                % already in dB
end
end

function x = nanZero(x)
% Zero (no data / no detection) -> NaN, so it is left out of dB plots.
x = double(x);
x(x <= 0) = NaN;
end
