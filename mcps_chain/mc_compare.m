function R = mc_compare(mine, rows, ref, refInfo, name, P, shift)
%MC_COMPARE Compare the output of one block with the matching log lane.
%
%   mc_compare()                                         print the table header
%   R = mc_compare(mine, rows, ref, refInfo, name, P)
%   R = mc_compare(mine, rows, ref, refInfo, name, P, shift)   use this range shift
%
%   mine    : my output, 2-D (pulses x cells) or 3-D (rows x cells x bins)
%   rows    : row bookkeeping of mine (output of the blocks, see mc_rows)
%   ref     : the log lane, refInfo its info (pulseSeq pairs the rows)
%
%   1. rows   pair the rows of the same pulse (pulseSeq); my rows with an
%             incomplete history are left out
%   2. range  shift that lines up the mean range profiles: ref(c) = mine(c + shift)
%   3. gain   least-squares gain |ref| = g |mine|
%   4. numbers on the magnitudes after the shift:
%        corr     correlation of the magnitudes
%        corr dB  correlation of the dB values
%        slope    slope of dB(log) against dB(mine): 1 = same law, 2 = the log is a
%                 power and mine an amplitude (0.5 = the reverse)
%        NMSE     error energy after the gain / log energy [dB]
%        coh      complex lanes: |sum(mine conj(log))| / (|mine| |log|),
%                 1 = equal up to one complex constant
%        det      CFAR lanes: detections in both / only mine / only the log
%                 (gain, corr ... then use only the cells detected by both)
%   3-D lanes are also compared after the max over the Doppler bins ([max]).

if nargin == 0
    fprintf('  %-32s %6s %8s %6s %9s %7s %8s %6s %9s %7s   %s\n', 'block', 'rows', 'paired', ...
            'shift', 'gain[dB]', 'corr', 'corr dB', 'slope', 'NMSE[dB]', 'coh', 'detections');
    return
end
R = struct('name', name, 'rows', 0, 'shift', NaN, 'gainDb', NaN, 'corr', NaN, ...
           'corrDb', NaN, 'slope', NaN, 'nmseDb', NaN, 'coh', NaN, 'det', [], 'max', []);
if isempty(mine) || isempty(ref)
    fprintf('  %-32s not available, skipped\n', name);
    return
end
C = P.compare;

% ---- 1. rows -----------------------------------------------------------------------------
[im, ir, how] = mc_match_rows(rows, refInfo, size(ref, 1));
if isempty(im)
    fprintf('  %-32s no common pulses (check pulseSeq / P.doppler.logPulseShift)\n', name);
    return
end
R.rows = numel(im);
if numel(im) > C.maxRows                                      % evenly spread subset
    pick = round(linspace(1, numel(im), C.maxRows));
    im = im(pick);
    ir = ir(pick);
end
nc = min(size(mine, 2), size(ref, 2));
nb = min(size(mine, 3), size(ref, 3));
A  = mine(im, 1:nc, 1:nb);
B  = ref(ir, 1:nc, 1:nb);
isDet = nnz(B) < 0.5 * numel(B);                              % mostly zeros = CFAR map

% ---- 2. range shift ----------------------------------------------------------------------
if nargin < 7 || isempty(shift)
    shift = mc_range_shift(profile(A), profile(B), 1:nc, C.maxShift);
end
A = shiftCells(A, shift);
R.shift = shift;
ok = (1:nc) + shift >= 1 & (1:nc) + shift <= nc;              % cells present in both

% ---- 3./4. gain and numbers ------------------------------------------------------------------
m = numbers(A(:, ok, :), B(:, ok, :), isDet);
R = copyNumbers(R, m);
printLine(name, R.rows, how, shift, m);
if nb > 1
    Am = max(abs(A), [], 3);
    Bm = max(abs(B), [], 3);
    R.max = numbers(Am(:, ok), Bm(:, ok), isDet);
    printLine([name ' [max]'], R.rows, how, shift, R.max);
else
    Am = A;
    Bm = B;
end

% ---- figure ----------------------------------------------------------------------------------
if C.plot
    g  = 10^(m.gainDb / 20);
    a  = 20*log10(double(abs(Am)) * g);
    b  = 20*log10(double(abs(Bm)));
    cl = mc_clim(b);
    figure('Name', name, 'NumberTitle', 'off');
    subplot(2, 2, 1);
    mc_show(Am * g, sprintf('mine x gain, shifted %+d cells', shift), cl);
    subplot(2, 2, 2);
    mc_show(Bm, 'log', cl);
    subplot(2, 2, 3);
    if isDet
        map = zeros(size(Am));                              % 0 none, 1 both, 2 mine, 3 log
        map(Am > 0 & Bm > 0)  = 1;
        map(Am > 0 & Bm == 0) = 2;
        map(Am == 0 & Bm > 0) = 3;
        [map, xc, yr] = mc_thin(map, 800, 600);          % differences drawn on top
        imagesc(xc, yr, map, [0 3]);
        axis xy;
        colormap(gca, [1 1 1; 0 0.6 0; 0.85 0 0; 0 0 0.85]);
        title('detections: green both, red only mine, blue only log');
    else
        d = a - b;
        d(~isfinite(d)) = 0;
        imagesc(d, [-10 10]);
        axis xy;
        colorbar;
        title('mine x gain - log [dB]');
    end
    xlabel('range cell');
    ylabel('row (pulse)');
    subplot(2, 2, 4);
    if isDet
        plot(1:nc, sum(Am > 0, 1), 1:nc, sum(Bm > 0, 1));
        ylabel('detections per cell');
    else
        plot(1:nc, 10*log10(profile(A) * g^2), 1:nc, 10*log10(profile(B)));
        ylabel('mean power [dB]');
    end
    grid on;
    legend('mine', 'log');
    xlabel('range cell');
    title(sprintf('range profile, %d rows paired by %s', R.rows, how));
end
end

% ==============================================================================================
function p = profile(X)
% mean power along range (over rows and bins), 1 x cells, double
p = double(squeeze(mean(mean(abs(X).^2, 1), 3)));
p = p(:).';
end

function Y = shiftCells(X, s)
% Y(:, c, :) = X(:, c + s, :), zeros where c + s is outside
nc = size(X, 2);
Y  = zeros(size(X), 'like', X);
if ~isreal(X)
    Y = complex(Y);
end
c  = (1:nc);
ok = c + s >= 1 & c + s <= nc;
Y(:, ok, :) = X(:, c(ok) + s, :);
end

function m = numbers(A, B, isDet)
a = abs(A(:));
b = abs(B(:));
m.det = [];
if isDet                                    % CFAR map: values only where both detect
    m.det = [nnz(a > 0 & b > 0), nnz(a > 0 & b == 0), nnz(a == 0 & b > 0)];
    both = a > 0 & b > 0;
    A = A(both);
    B = B(both);
    a = a(both);
    b = b(both);
end
m.gainDb = 20*log10(sum(a .* b, 'double') / sum(a.^2, 'double'));
g        = 10^(m.gainDb / 20);
m.corr   = corrOf(a, b);
both     = a > 0 & b > 0;
[m.corrDb, m.slope] = corrOf(20*log10(a(both)), 20*log10(b(both)));
m.nmseDb = 10*log10(sum((g * a - b).^2, 'double') / sum(b.^2, 'double'));
m.coh    = NaN;
if ~isreal(A) && ~isreal(B)
    m.coh = abs(sum(A(:) .* conj(B(:)), 'double')) / ...
            sqrt(sum(a.^2, 'double') * sum(b.^2, 'double'));
end
end

function [r, slope] = corrOf(x, y)
% correlation and least-squares slope of y against x (double sums)
r = NaN;
slope = NaN;
if numel(x) < 2
    return
end
x = x - sum(x, 'double') / numel(x);
y = y - sum(y, 'double') / numel(y);
sxy = sum(x .* y, 'double');
sxx = sum(x .* x, 'double');
syy = sum(y .* y, 'double');
r = sxy / sqrt(sxx * syy);
slope = sxy / sxx;
end

function R = copyNumbers(R, m)
f = fieldnames(m);
for i = 1:numel(f)
    R.(f{i}) = m.(f{i});
end
end

function printLine(name, nRows, how, shift, m)
fprintf('  %-32s %6d %8s %+6d %+9.2f %7.4f %8.4f %6.2f %+9.1f %7.4f', name, nRows, how, ...
        shift, m.gainDb, m.corr, m.corrDb, m.slope, m.nmseDb, m.coh);
if ~isempty(m.det)
    fprintf('   both %d / only mine %d / only log %d', m.det);
end
fprintf('\n');
end
