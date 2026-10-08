function [X, demo] = if_demo_data()
%IF_DEMO_DATA Synthetic raw data to try main_if.m.
%
%   [X, demo] = if_demo_data()
%
%   X     : 400 PRIs x 1500 range cells, complex, noise power 1
%   demo  : .short, .long  the two transmitted codes: a short Barker 13 and a long
%                          64-chip code; both are sent in every PRI
%           .ifDb          true improvement factor of the simulated clutter
%           .testIfDb      true improvement factor of the test target
%
%   Ground clutter is strong at near range and falls with range. It moves a
%   little (Gaussian Doppler spectrum, width sigma = 2 % of the PRF), so the
%   canceler cannot remove all of it. One moving target sits inside the clutter.
%   A test target (as from a target generator) sits at cell 1000 in every PRI,
%   a ring on the PPI. It does not move; only a small phase jitter limits its IF.

rng(1);
nPri  = 400;
nCell = 1500;
sigma = 0.02;                                            % clutter Doppler spread / PRF

short = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];                % Barker 13
long  = exp(2j*pi*rand(1, 64));                          % 64-chip random-phase code
tx    = [short, zeros(1, 30), long];                     % both pulses in every PRI

% One clutter scatterer per range cell, 50 dB at cell 1 and 0.1 dB weaker per cell.
% Each one fluctuates from PRI to PRI with a Gaussian Doppler spectrum.
f     = mod((0:nPri-1)' / nPri + 0.5, 1) - 0.5;          % Doppler / PRF, FFT order
A     = ifft(exp(-f.^2 / (4*sigma^2)) .* fft(randn(nPri, nCell) + 1j*randn(nPri, nCell)));
A     = A / sqrt(mean(abs(A(:)).^2));                    % unit power
scat  = A .* 10.^((50 - 0.1*(0:nCell-1)) / 20);

m = (0:nPri-1)';
scat(:, 300) = scat(:, 300) + 100 * exp(2j*pi*0.3*m);  % moving target: 40 dB, Doppler 0.3 PRF
jitter = 10^(-70/20);                                    % rms phase jitter [rad]: IF = 1/jitter^2
scat(:, 1000) = scat(:, 1000) + 1e4 * exp(1j * jitter * randn(nPri, 1));   % test target, 80 dB

X = filter(tx, 1, scat, [], 2);                          % raw echo of both pulses
X = X + (randn(nPri, nCell) + 1j*randn(nPri, nCell)) / sqrt(2);   % receiver noise, power 1

rho           = @(k) exp(-2 * (pi * sigma * k)^2);       % clutter correlation at lag k
demo.short    = short;
demo.long     = long;
demo.ifDb     = 10*log10(6 / (6 - 8*rho(1) + 2*rho(2)));
demo.testIfDb = -20*log10(jitter);
end
