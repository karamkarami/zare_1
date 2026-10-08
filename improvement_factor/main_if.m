%% MAIN_IF  Improvement factor of a 3-pulse MTI canceler, measured on raw complex data
%
%   Data  : complex matrix, one row per PRI, one column per range cell.
%   Steps : load -> noise level -> find the clutter cells -> canceler -> IF -> plots
%
%   The IF is only as good as the choice of clutter cells. The script takes every
%   strong echo near zero Doppler: ground clutter, but also a test target (a ring on
%   the PPI), a delay line or transmitter leakage. Check the table of clutter
%   segments it prints, and choose the cells with clutterCells / skipCells.
%   A steady test target is a good "clutter" for the IF of the radar itself: it
%   does not move, so all that the canceler leaves of it comes from the radar.
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
canc = @(Z) w(1) * Z(1+2*s:end, :) + w(2) * Z(1+s:end-s, :) + w(3) * Z(1:end-2*s, :);
Y    = canc(X);
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

% Cells far more stable than the rest: a test target, a delay line or leakage
stable = c & 10*log10(G * Pin ./ Pout) > median(ifCellDb(c & ~isnan(ifCellDb))) + outlierDb;
if any(stable)
    fprintf('Note: %d cells between cells %d and %d are far more stable than the rest\n', ...
            sum(stable), min(cells(stable)), max(cells(stable)));
    fprintf('      (test target, delay line, transmitter leakage?).\n');
end

% IF of a group of cells m (true/false mask or indices) from an output power P
% whose noise part is Pn in every cell: with the noise removed (ifOf), or without
% removing it (lowOf, a lower bound). okOf: the clutter left is above 10 % of the noise.
nOf   = @(m) sum(m > 0);
ifOf  = @(m, P, Pn) 10*log10(G * sum(Cin(m)) / max(sum(P(m)) - Pn * nOf(m), eps));
lowOf = @(m, P)     10*log10(G * sum(Pin(m)) / sum(P(m)));
okOf  = @(m, P, Pn) sum(P(m)) - Pn * nOf(m) > 0.1 * Pn * nOf(m);

% Steady echo of every cell: the same echo in every PRI, turning with a constant
% Doppler (phase step angle(R1) from one PRI to the next). It is the mean over the
% PRIs after the Doppler is taken out. For a test target nearly all of the power
% is steady, for moving ground clutter only a little.
ref = zeros(size(X));
for q = 1:s                                              % with priStep = 2: each pulse alone
    rows = q:s:nPri;
    rot  = exp(1j * (0:numel(rows)-1)' * angle(R1));
    ref(rows, :) = mean(X(rows, :) .* conj(rot), 1) .* rot;
end
steadyOf = @(m) sum(sum(abs(ref(:, m)).^2)) / (nPri * sum(Pin(m)));

% What the canceler leaves of a steady echo, in three parts:
%   Doppler   : the canceler output of the steady echo itself (its constant Doppler)
%   amplitude : the changes from pulse to pulse in phase with the echo
%   phase     : the changes at 90 degrees to the echo (phase noise, timing jitter)
% The last two each carry half of the output noise.
Yd    = canc(X - ref);                                   % canceler output of the changes
u     = ref(1+s:end-s, :) ./ abs(ref(1+s:end-s, :));     % phase of the steady echo
PoutD = mean(abs(canc(ref)).^2, 1);                      % Doppler part
PoutA = mean(real(Yd .* conj(u)).^2, 1);                 % amplitude part
PoutP = mean(imag(Yd .* conj(u)).^2, 1);                 % phase part

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
shown = order(1:min(8, nSeg))';
fprintf('\nClutter segments, strongest first:\n');
fprintf('     cells        max CNR   power share   steady       IF\n');
for k = shown
    m = seg{k};
    if okOf(m, Pout, G * noise)
        ifTxt = sprintf('   %5.1f dB', ifOf(m, Pout, G * noise));
    else
        ifTxt = sprintf('>= %5.1f dB', lowOf(m, Pout));
    end
    fprintf('  %5d - %-5d  %5.1f dB  %8.1f %%   %5.1f %%   %s\n', cells(m(1)), cells(m(end)), ...
            max(cnrDb(m)), 100 * share(k), 100 * steadyOf(m), ifTxt);
end
if nSeg > 8
    fprintf('  (+ %d smaller segments)\n', nSeg - 8);
end

% Steady segments (test target, delay line): the IF that each part alone would
% give. 1/IF = 1/IFdoppler + 1/IFamplitude + 1/IFphase. A test target that is not
% locked to the radar has a small Doppler f, which alone limits the IF to
% 6 / (16 sin^4(pi f)).
shownSteady = shown(cellfun(steadyOf, seg(shown)) > 0.9);
if ~isempty(shownSteady)
    fprintf('\nSteady segments (more than 90 %% of the power is the same echo in every PRI):\n');
    fprintf('     cells       Doppler/PRF    Doppler part  amplitude part    phase part\n');
    for k = shownSteady
        m   = seg{k};
        txt = sprintf('     %6.1f dB', min(10*log10(G * sum(Cin(m)) / sum(PoutD(m))), 999.9));
        for P = {PoutA, PoutP}
            if okOf(m, P{1}, G * noise / 2)
                txt = [txt, sprintf('      %5.1f dB', ifOf(m, P{1}, G * noise / 2))];
            else
                txt = [txt, sprintf('   >= %5.1f dB', lowOf(m, P{1}))];
            end
        end
        fprintf('  %5d - %-5d   %+9.6f %s\n', cells(m(1)), cells(m(end)), angle(sum(R1(m))) / (2*pi), txt);
    end
end

fprintf('\n');
if okOf(c, Pout, G * noise)
    ifDb = ifOf(c, Pout, G * noise);
    caDb = ifDb - 10*log10(G);
    fprintf('Clutter attenuation  CA = %.1f dB\n', caDb);
    fprintf('Improvement factor   IF = CA + 10log10(%d) = %.1f dB   (%d clutter cells)\n', G, ifDb, nCl);
    fprintf('(without noise removal IF = %.1f dB, a lower bound)\n', lowOf(c, Pout));
else
    ifDb = lowOf(c, Pout);
    fprintf('Improvement factor   IF >= %.1f dB\n', ifDb);
    fprintf('The clutter left after the canceler is below the noise, so only a lower bound\n');
    fprintf('can be measured. Use stronger clutter cells (raise minCnrDb).\n');
end
ifMedDb = median(ifCellDb(c & ~isnan(ifCellDb)));
fprintf('IF of single clutter cells: median %.1f dB\n', ifMedDb);
if ifDb > ifMedDb + 6
    k = order(1);
    fprintf('\nThe total IF is far above the typical cell: it is set by the strongest segment,\n');
    fprintf('cells %d-%d (%.0f %% of the clutter power). To measure on ground clutter, leave\n', ...
            cells(seg{k}(1)), cells(seg{k}(end)), 100 * share(k));
    fprintf('it out with skipCells; to measure on that echo alone, set clutterCells to it.\n');
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
