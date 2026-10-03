function r = rsp_compare(ref, test, name, varargin)
%RSP_COMPARE Compare a block output with a reference lane (e.g. the MCPS log).
%
%   r = rsp_compare(ref, test, name, 'option', value, ...)
%
%   ref, test : arrays of the same layout (pulses x range [x doppler]);
%               sizes may differ, the common part is compared
%   name      : label printed in the summary line
%
%   Options
%     'align'     true      search the pulse / range shift of test vs ref
%     'maxLag'    [8 64]    largest pulse and range shift searched
%     'tolerance' 1         detections: range cells accepted around a hit
%     'plot'      false     show ref, scaled test and their difference
%     'quiet'     false     do not print
%
%   Output r
%     .lag        [pulse range] shift: test(i + lag) matches ref(i)
%     .size       compared size
%     .gain       least-squares gain g minimising |ref - g*test|
%     .gainDb     20*log10(|g|)
%     .corr       |<ref, test>| / (|ref| |test|), 1 = identical up to a gain
%     .nmseDb     |ref - g*test|^2 / |ref|^2 in dB (-Inf = identical)
%   For detection maps (logical, or a CFAR lane that is 0 off detections):
%     .refCount .testCount .matched .missed .extra .jaccard

opt = struct('align', true, 'maxLag', [8 64], 'tolerance', 1, 'plot', false, 'quiet', false);
for i = 1:2:numel(varargin)
    opt.(varargin{i}) = varargin{i+1};
end
if nargin < 3 || isempty(name)
    name = 'lane';
end

r = struct('name', name, 'lag', [0 0], 'size', [], 'gain', NaN, 'gainDb', NaN, ...
           'corr', NaN, 'nmseDb', NaN);
if isempty(ref) || isempty(test)
    if ~opt.quiet
        fprintf('%-28s  skipped (empty input)\n', name);
    end
    return
end

isDet = islogical(ref) || islogical(test);

% --- shift between the two arrays -------------------------------------------------
if opt.align
    for dim = 1:2
        r.lag(dim) = bestLag(profile(ref, dim), profile(test, dim), opt.maxLag(dim));
    end
end
[ref, test] = overlap(ref, test, r.lag);
r.size = size(ref);

if isDet || isDetectionMap(ref, test)
    % --- detection maps ---------------------------------------------------------------
    a = ref ~= 0;
    b = test ~= 0;
    da = dilateRange(a, opt.tolerance);
    db = dilateRange(b, opt.tolerance);
    r.refCount  = nnz(a);
    r.testCount = nnz(b);
    r.matched   = nnz(a & db);
    r.missed    = nnz(a & ~db);
    r.extra     = nnz(b & ~da);
    r.jaccard   = nnz(a & b) / max(nnz(a | b), 1);
    if ~opt.quiet
        fprintf(['%-28s  lag [%d %d]  size %-16s  ref %7d  test %7d  matched %7d  ' ...
                 'missed %7d  extra %7d  Jaccard %.3f\n'], name, r.lag, mat2str(r.size), ...
                r.refCount, r.testCount, r.matched, r.missed, r.extra, r.jaccard);
    end
end

if ~islogical(ref) && ~islogical(test)
    % --- sample values ----------------------------------------------------------------
    % one page (2-D slice) at a time: no full-size double copies
    nPg = numel(ref) / (size(ref, 1) * size(ref, 2));
    xx = 0;  yy = 0;  xy = 0;
    for k = 1:nPg
        x  = double(reshape(ref(:, :, k), [], 1));
        y  = double(reshape(test(:, :, k), [], 1));
        xx = xx + real(x' * x);
        yy = yy + real(y' * y);
        xy = xy + y' * x;
    end
    r.gain = xy / max(yy, realmin);
    ee = 0;
    for k = 1:nPg
        x  = double(reshape(ref(:, :, k), [], 1));
        y  = double(reshape(test(:, :, k), [], 1));
        ee = ee + sum(abs(x - r.gain*y).^2);
    end
    r.gainDb = 20*log10(abs(r.gain));
    r.corr   = abs(xy) / sqrt(max(xx*yy, realmin));
    r.nmseDb = 10*log10(ee / max(xx, realmin));
    if ~opt.quiet
        fprintf('%-28s  lag [%d %d]  size %-16s  gain %8.2f dB  corr %.6f  NMSE %8.2f dB\n', ...
                name, r.lag, mat2str(r.size), r.gainDb, r.corr, r.nmseDb);
    end
end

if opt.plot
    plotCompare(ref, test, r);
end
end

% ======================================================================================
function p = profile(x, dim)
% Power of x (dim = 1: per pulse, dim = 2: per range cell), summed over
% every other dimension; one page at a time.
nPg = numel(x) / (size(x, 1) * size(x, 2));
p   = zeros(size(x, dim), 1);
for k = 1:nPg
    a = abs(double(x(:, :, k))).^2;
    p = p + reshape(sum(a, 3 - dim), [], 1);
end
end

function lag = bestLag(a, b, maxLag)
% Shift of b against a with the largest correlation of the centred profiles.
a = a - mean(a);
b = b - mean(b);
n = 2^nextpow2(numel(a) + numel(b));
c = ifft(conj(fft(a, n)) .* fft(b, n));            % c(l) = sum a(i) b(i + l)
l = -maxLag:maxLag;
v = real(c(mod(l, n) + 1));
[~, k] = max(v);
lag = l(k);
end

function [a, b] = overlap(a, b, lag)
% Common part of a and b after shifting b by lag (pulse, range).
sa = size(a);
sb = size(b);
nd = max(numel(sa), numel(sb));
sa(end+1:nd) = 1;
sb(end+1:nd) = 1;
ia = cell(1, nd);
ib = cell(1, nd);
for d = 1:nd
    l = 0;
    if d <= 2
        l = lag(d);
    end
    n = min(sa(d) - max(-l, 0), sb(d) - max(l, 0));
    ia{d} = max(-l, 0) + (1:n);
    ib{d} = max(l, 0) + (1:n);
end
a = a(ia{:});
b = b(ib{:});
end

function tf = isDetectionMap(a, b)
% A CFAR lane stored as values on detections and 0 elsewhere.
tf = isreal(a) && isreal(b) && nnz(a) < 0.2*numel(a) && nnz(b) < 0.2*numel(b);
end

function d = dilateRange(x, t)
% Mark cells within t range cells (dimension 2) of a true cell.
d = x;
for s = 1:t
    d(:, 1+s:end, :) = d(:, 1+s:end, :) | x(:, 1:end-s, :);
    d(:, 1:end-s, :) = d(:, 1:end-s, :) | x(:, 1+s:end, :);
end
end

function plotCompare(ref, test, r)
% Reference, gain-matched test and difference (max over Doppler for 3-D lanes).
if islogical(ref) || islogical(test)
    g = 1;
else
    g = r.gain;
end
toDb = @(x) 20*log10(max(abs(double(x)), realmin));
A = max(toDb(ref), [], 3);
B = max(toDb(g * test), [], 3);
E = max(toDb(double(ref) - g*double(test)), [], 3);
lim = rsp_db_limits(A);
figure('Name', ['Compare: ' r.name], 'Color', 'w');
subplot(1, 3, 1); imagesc(A, lim); axis xy; colorbar; title('reference [dB]');
xlabel('range cell'); ylabel('pulse');
subplot(1, 3, 2); imagesc(B, lim); axis xy; colorbar; title('test x gain [dB]');
xlabel('range cell');
subplot(1, 3, 3); imagesc(E, lim); axis xy; colorbar;
title(sprintf('difference [dB], NMSE %.1f dB', r.nmseDb));
xlabel('range cell');
colormap(jet);
end
