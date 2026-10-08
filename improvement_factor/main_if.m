%% MAIN_IF  Improvement factor of a 3-pulse MTI canceler, measured on raw complex data
%
%   Data  : complex matrix, one row per PRI, one column per range cell.
%   Steps : load -> noise level -> find the clutter cells -> canceler -> IF -> plots
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
% Clipped samples pile up at the largest value.
iq   = abs([real(X(:)); imag(X(:))]);
nTop = sum(iq >= 0.999 * max(iq));
if nTop > 10
    warning('%d samples sit at the largest value: the ADC is probably saturated, the IF will be limited.', nTop);
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
strong = cnrDb >= minCnrDb;
if ~any(strong)
    error('No cell is %g dB above the noise. Lower minCnrDb, or check firstCell and the data.', minCnrDb);
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

% Cells far more stable than the rest. At the very start of the range this is
% transmitter leakage, which makes the IF look too good: move firstCell past it.
stable = c & 10*log10(G * Pin ./ Pout) > median(ifCellDb(c & ~isnan(ifCellDb))) + outlierDb;
if any(stable)
    fprintf('Note: %d cells between cells %d and %d are far more stable than the rest.\n', ...
            sum(stable), min(cells(stable)), max(cells(stable)));
    fprintf('      If they are at the start of the range (transmitter leakage), set firstCell after them.\n');
end

ifLowDb = 10*log10(G * sum(Pin(c)) / sum(Pout(c)));      % noise not removed: a lower bound
residue = sum(Cout(c)) / (G * noise * nCl);              % clutter left / noise, at the output

fprintf('\n');
if residue > 0.1
    caDb = 10*log10(sum(Cin(c)) / sum(Cout(c)));
    ifDb = caDb + 10*log10(G);
    fprintf('Clutter attenuation  CA = %.1f dB\n', caDb);
    fprintf('Improvement factor   IF = CA + 10log10(%d) = %.1f dB   (%d clutter cells)\n', G, ifDb, nCl);
    fprintf('(without noise removal IF = %.1f dB, a lower bound)\n', ifLowDb);
else
    ifDb = ifLowDb;
    fprintf('Improvement factor   IF >= %.1f dB\n', ifLowDb);
    fprintf('The clutter left after the canceler is below the noise, so only a lower bound\n');
    fprintf('can be measured. Use stronger clutter cells (raise minCnrDb).\n');
end
fprintf('IF of single clutter cells: median %.1f dB\n', median(ifCellDb(c & ~isnan(ifCellDb))));
if exist('demo', 'var')
    fprintf('Demo: true IF of the simulated clutter = %.1f dB\n', demo.ifDb);
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
