%% MAIN_IF  Improvement factor of a 3-pulse MTI canceler, measured on raw complex data
%
%   Data  : complex matrix, one row per PRI, one column per range cell.
%   Steps : load -> noise level -> find the clutter cells -> canceler -> IF -> plots
%
%   The IF is only as good as the choice of clutter cells. The script takes every
%   strong echo near zero Doppler, so a test target (a ring on the PPI), a delay
%   line or transmitter leakage looks like clutter to it. Check the table of
%   clutter segments it prints, and set clutterCells / skipCells when needed.
%
%   Canceler            y(m) = x(m) - 2 x(m-1) + x(m-2)
%   Improvement factor  IF = (S/C)out / (S/C)in = G * Cin / Cout
%                       G = 1^2 + 2^2 + 1^2 = 6 is the canceler gain for noise and
%                       the mean gain for targets of all Doppler frequencies.
%
%   A decoder is not needed: it filters along range, the canceler along the PRIs,
%   and both are linear, so their order does not change the canceler output.
%   Set decoderTaps to check this on your own data.

clear; clc; close all;

%% 1. Settings ----------------------------------------------------------------------
dataFile     = '';     % .mat file with your data; '' = synthetic demo (if_demo_data.m)
dataVar      = '';     % variable in the file; '' = the largest numeric variable
pulsesInRows = true;   % true: rows = PRIs, columns = range cells; false: the opposite
pris         = [];     % PRIs to use, e.g. 1:300; [] = all
firstCell    = 1;      % first range cell to use: skip the transmitter leakage / blanking
clutterCells = [];     % where the clutter is, one [first last] per row, e.g. [1 400]; [] = anywhere
skipCells    = [];     % cells that are not clutter, e.g. a test target: [990 1110]
decoderTaps  = [];     % optional decoder taps along range, e.g. conj(fliplr(code)); [] = raw data
% demo check: [~, dd] = if_demo_data(); decoderTaps = conj(fliplr(dd.long));
priStep      = 1;      % 1 = cancel consecutive PRIs; 2 = every other PRI
                       %     (only if the two pulses alternate from PRI to PRI)
minCnrDb     = 20;     % a clutter cell is at least this many dB above the noise ...
maxDoppler   = 0.05;   % ... and its mean Doppler is below this fraction of the PRF
outlierDb    = 10;     % drop clutter cells whose IF is this far below the median (moving objects)

%% 2. Load the data -----------------------------------------------------------------
if isempty(dataFile)
    [X, demo] = if_demo_data();
    fprintf('Synthetic demo data (if_demo_data.m)\n');
else
    S = load(dataFile);
    if isempty(dataVar)                                  % largest numeric variable
        names = fieldnames(S);
        sizes = cellfun(@(f) isnumeric(S.(f)) * numel(S.(f)), names);
        [~, k] = max(sizes);
        dataVar = names{k};
    end
    X = S.(dataVar);
    fprintf('File %s, variable "%s"\n', dataFile, dataVar);
end
if ~pulsesInRows
    X = X.';                                             % .' = transpose without conjugate
end
X = double(X);
if ~isempty(pris)
    X = X(pris, :);
end

% ADC saturation makes the system nonlinear before any digital processing.
% Clipped samples pile up at exactly the same largest (or smallest) value.
iq   = [real(X(:)); imag(X(:))];
nTop = max(sum(iq == max(iq)), sum(iq == min(iq)));
if nTop > 10
    warning('%d samples sit exactly at the largest value: is the ADC saturated? Then the IF is limited.', nTop);
end

X     = X(:, firstCell:end);
cells = firstCell : firstCell + size(X, 2) - 1;          % range cell numbers
if ~isempty(decoderTaps)
    h = decoderTaps(:).';
    X = filter(h, 1, X, [], 2) / norm(h);                % same taps on every PRI, noise level kept
end
[nPri, nCell] = size(X);
fprintf('%d PRIs x %d range cells (cells %d..%d)\n', nPri, nCell, cells(1), cells(end));

%% 3. Noise level -------------------------------------------------------------------
% Clutter covers only part of the range, so the quietest cells hold noise only:
% quiet = cells less than 1.8 dB (x 1.5) above the 10th percentile of the power profile.
Pin   = mean(abs(X).^2, 1);                              % mean power of every range cell
Ps    = sort(Pin(Pin > 0));
p10   = Ps(max(1, round(0.1 * numel(Ps))));
quiet = Pin > 0 & Pin < 1.5 * p10;

% A receiver DC offset looks like perfectly stable clutter in every cell: remove it.
Q     = X(:, quiet);
dc    = mean(Q(:));
X     = X - dc;
Pin   = mean(abs(X).^2, 1);
noise = median(Pin(quiet));
fprintf('Noise power %.4g (%.1f dB), DC offset %.1f dB relative to the noise (removed)\n', ...
        noise, 10*log10(noise), 10*log10(abs(dc)^2 / noise));

%% 4. Clutter cells -----------------------------------------------------------------
cnrDb  = 10*log10(max(Pin / noise - 1, eps));            % clutter-to-noise ratio of every cell

area = true(1, nCell);                                   % where clutter may be
if ~isempty(clutterCells)
    area(:) = false;
    for k = 1:size(clutterCells, 1)
        area(cells >= clutterCells(k, 1) & cells <= clutterCells(k, 2)) = true;
    end
end
for k = 1:size(skipCells, 1)
    area(cells >= skipCells(k, 1) & cells <= skipCells(k, 2)) = false;
end

strong = cnrDb >= minCnrDb & area;
if ~any(strong)
    error('No cell is %g dB above the noise. Lower minCnrDb, or check firstCell, clutterCells and the data.', minCnrDb);
end

% Correlation of the strong cells from one PRI to the next (lag 1) and two PRIs on
% (lag 2). Stable clutter gives values close to 1 at both lags. A low lag 1 with a
% high lag 2 means that the PRIs alternate between two different pulses.
lagCorr = @(Z, k) abs(mean(mean(Z(1+k:end, :) .* conj(Z(1:end-k, :))))) / mean(mean(abs(Z).^2));
rho1 = lagCorr(X(:, strong), 1);
rho2 = lagCorr(X(:, strong), 2);
fprintf('PRI-to-PRI correlation of the strong cells: lag 1 = %.4f, lag 2 = %.4f\n', rho1, rho2);
if priStep == 1 && rho2 > rho1 + 0.1
    warning('Lag 2 is much more correlated than lag 1: the pulses seem to alternate. Set priStep = 2.');
end

% Mean Doppler of every cell (pulse pair): ground clutter sits near zero Doppler,
% moving things (aircraft, cars, birds, rain) do not.
s   = priStep;
R1  = mean(X(1+s:end, :) .* conj(X(1:end-s, :)), 1);
dop = angle(R1) / (2*pi);                                % mean Doppler / PRF seen by the canceler
c   = strong & abs(dop) <= maxDoppler;                   % the clutter cells
nCl = sum(c);
if nCl == 0
    error('No strong cell is near zero Doppler. Check priStep, or raise maxDoppler.');
end
fprintf('Clutter cells: %d between cells %d and %d, median CNR %.1f dB (%d strong cells rejected as moving)\n', ...
        nCl, min(cells(c)), max(cells(c)), median(cnrDb(c)), sum(strong & ~c));

%% 5. 3-pulse canceler --------------------------------------------------------------
w    = [1 -2 1];
Y    = w(1) * X(1+2*s:end, :) + w(2) * X(1+s:end-s, :) + w(3) * X(1:end-2*s, :);
G    = sum(abs(w).^2);                                   % = 6
Pout = mean(abs(Y).^2, 1);

%% 6. Improvement factor ------------------------------------------------------------
Cin  = Pin  - noise;                                     % clutter power at the input
Cout = Pout - G * noise;                                 % clutter left at the output

% IF of every single clutter cell (NaN where the clutter left is below the noise)
r = G * Cin ./ Cout;
r(~c | Cout <= 0) = NaN;
ifCellDb = 10*log10(r);

% A moving object (aircraft, car, bird) inside a clutter cell passes the canceler
% and makes the cell look like very unstable clutter. Drop such cells.
odd = ifCellDb < median(ifCellDb(~isnan(ifCellDb))) - outlierDb;
c   = c & ~odd;
nCl = sum(c);
fprintf('%d clutter cells dropped: IF more than %g dB below the median (moving objects)\n', sum(odd), outlierDb);

% Cells far more stable than the rest: transmitter leakage, a test target or a
% delay line. They are not clutter and make the IF look too good.
stable = c & 10*log10(G * Pin ./ Pout) > median(ifCellDb(c & ~isnan(ifCellDb))) + outlierDb;
if any(stable)
    fprintf('Note: %d cells between cells %d and %d are far more stable than the rest.\n', ...
            sum(stable), min(cells(stable)), max(cells(stable)));
    fprintf('      Transmitter leakage, a test target, a delay line? If so, leave them out (skipCells).\n');
end

% IF of a group of cells m (true/false mask or indices)
ifOf  = @(m) 10*log10(G * sum(Cin(m)) / sum(Cout(m)));  % noise removed
lowOf = @(m) 10*log10(G * sum(Pin(m)) / sum(Pout(m)));  % noise not removed: a lower bound
okOf  = @(m) sum(Cout(m)) > 0.1 * G * noise * sum(m > 0);   % clutter left is above 10 % of the noise

% Clutter segments: neighbouring clutter cells (gaps up to 3 cells). The total IF
% is weighted by power, so the strongest segment sets it.
idx   = find(c);
brk   = [0, find(diff(idx) > 3), numel(idx)];
nSeg  = numel(brk) - 1;
seg   = cell(nSeg, 1);
share = zeros(nSeg, 1);
for k = 1:nSeg
    seg{k}   = idx(brk(k)+1 : brk(k+1));
    share(k) = sum(Cin(seg{k})) / sum(Cin(c));
end
[~, order] = sort(share, 'descend');
fprintf('\nClutter segments, strongest first:\n');
fprintf('     cells        max CNR   power share      IF\n');
for k = order(1:min(8, nSeg))'
    m = seg{k};
    if okOf(m)
        ifTxt = sprintf('   %5.1f dB', ifOf(m));
    else
        ifTxt = sprintf('>= %5.1f dB', lowOf(m));
    end
    fprintf('  %5d - %-5d  %5.1f dB  %8.1f %%   %s\n', cells(m(1)), cells(m(end)), max(cnrDb(m)), 100 * share(k), ifTxt);
end
if nSeg > 8
    fprintf('  (+ %d smaller segments)\n', nSeg - 8);
end

fprintf('\n');
if okOf(c)
    caDb = ifOf(c) - 10*log10(G);
    ifDb = ifOf(c);
    fprintf('Clutter attenuation  CA = %.1f dB\n', caDb);
    fprintf('Improvement factor   IF = CA + 10log10(%d) = %.1f dB   (%d clutter cells)\n', G, ifDb, nCl);
    fprintf('(without noise removal IF = %.1f dB, a lower bound)\n', lowOf(c));
else
    ifDb = lowOf(c);
    fprintf('Improvement factor   IF >= %.1f dB\n', ifDb);
    fprintf('The clutter left after the canceler is below the noise, so only a lower bound\n');
    fprintf('can be measured. Use stronger clutter cells (raise minCnrDb).\n');
end
ifMedDb = median(ifCellDb(c & ~isnan(ifCellDb)));
fprintf('IF of single clutter cells: median %.1f dB\n', ifMedDb);
if ifDb > ifMedDb + 6
    k = order(1);
    fprintf('\nThe total IF is far above the typical cell: it is set by the strongest segment,\n');
    fprintf('cells %d-%d (%.0f %% of the clutter power). If that is not ground clutter\n', ...
            cells(seg{k}(1)), cells(seg{k}(end)), 100 * share(k));
    fprintf('(test target, delay line, leakage), leave it out with skipCells.\n');
end
if exist('demo', 'var')
    fprintf('Demo: true IF of the clutter %.1f dB, of the test target (cell 1000) %.1f dB\n', demo.ifDb, demo.testIfDb);
end

%% 7. Plots -------------------------------------------------------------------------
figure('Name', 'Improvement factor', 'Color', 'w');

subplot(2, 2, 1);
imagesc(cells, 1:nPri, 10*log10(abs(X).^2 / noise));
caxis([-5, max(cnrDb) + 5]);
colorbar;
xlabel('range cell'); ylabel('PRI');
title('data [dB above noise]');

subplot(2, 2, 2);
plot(cells, 10*log10(Pin / noise), 'b', cells, 10*log10(Pout / noise), 'r');
hold on;
plot(cells(c), 10*log10(Pin(c) / noise), 'k.', 'MarkerSize', 6);
plot(cells([1 end]), [0 0], 'b:', cells([1 end]), 10*log10(G) * [1 1], 'r:');
grid on;
xlabel('range cell'); ylabel('dB above noise');
legend('canceler input', 'canceler output', 'clutter cells', 'input noise', 'output noise (+7.8 dB)', ...
       'Location', 'best');
xlim(cells([1 end]));
title('mean power of every range cell');

subplot(2, 2, 3);
kept    = ifCellDb;  kept(~c)      = NaN;
dropped = ifCellDb;  dropped(~odd) = NaN;
plot(cells, kept, 'k.', cells, dropped, 'mx');
hold on;
plot(cells([1 end]), ifDb * [1 1], 'r', 'LineWidth', 1.5);
grid on;
xlabel('range cell'); ylabel('IF [dB]');
legend('single cells', 'dropped cells', 'all clutter cells', 'Location', 'best');
xlim(cells([1 end]));
title(sprintf('improvement factor %.1f dB', ifDb));

% Doppler spectrum of the clutter cells before and after the canceler
subplot(2, 2, 4);
hannw = @(n) 0.5 - 0.5 * cos(2*pi*(0:n-1)' / n);
spec  = @(Z) fftshift(mean(abs(fft(Z .* hannw(size(Z, 1)), nPri)).^2, 2)) / sum(hannw(size(Z, 1)).^2);
f     = ((0:nPri-1)' - floor(nPri/2)) / nPri;            % Doppler / PRF seen by the canceler
gain  = 10*log10(16 * sin(pi * f).^4);                   % |1 - 2 z^-1 + z^-2|^2
Sin   = 10*log10(spec(X(1:s:end, c)) / noise);
plot(f, Sin, 'b', f, 10*log10(spec(Y(1:s:end, c)) / noise), 'r', f, gain, 'k--');
grid on;
xlim([-0.5 0.5]);
ylim([-40, max(Sin) + 10]);
xlabel('Doppler / PRF'); ylabel('dB');
legend('input [dB above noise]', 'output [dB above noise]', 'canceler gain [dB]', 'Location', 'best');
title('mean Doppler spectrum of the clutter cells');
