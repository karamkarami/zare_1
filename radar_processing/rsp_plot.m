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
%     'bin'      []      Doppler bin for the range cut ([] = strongest detection)
%     'range'    []      range cells shown, [first last] ([] = all)
%     'colormap' 'jet'

opt = struct('frame', [], 'bin', [], 'range', [], 'colormap', 'jet');
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
image2(panel + 3, max(pwDb(out.integral(:, rc, :)), [], 3), rc, ...
       sprintf('integral n=%d, max bin', out.P.nci.nFrames));

% CFAR detections: range x row, coloured by Doppler bin
subplot(2, 3, panel + 4);
L  = out.cfar.list;
in = L.range >= rc(1) & L.range <= rc(end);
scatter(L.range(in), L.frame(in), 6, L.bin(in), 'filled');
axis([rc(1) rc(end) 0.5 size(det, 1) + 0.5]);
cb = colorbar; ylabel(cb, 'Doppler bin');
xlabel('range cell'); ylabel('frame');
title(sprintf('%s-CFAR: %d detections', upper(out.P.cfar.type), numel(L.range)));

% range cut: integrated signal against the CFAR threshold
subplot(2, 3, panel + 5);
x = squeeze(out.integral(opt.frame, rc, opt.bin));
t = squeeze(out.cfar.threshold(opt.frame, rc, opt.bin));
h = squeeze(det(opt.frame, rc, opt.bin));
plot(rc, pwDb(x), 'Color', [0.2 0.4 0.8]); hold on
plot(rc, pwDb(t), 'r', 'LineWidth', 1);
plot(rc(h), pwDb(x(h)), 'ko', 'MarkerFaceColor', 'y');
hold off; grid on; xlim([rc(1) rc(end)]);
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
    case 'square', f = @(x) 10*log10(max(double(x), realmin));
    case 'linear', f = @(x) 20*log10(max(double(x), realmin));
    otherwise,     f = @(x) double(x);                % already in dB
end
end
