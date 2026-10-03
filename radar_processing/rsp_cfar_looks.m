function [nInt, nRef] = rsp_cfar_looks(P, replica)
%RSP_CFAR_LOOKS Effective looks seen by the CFAR, per Doppler bin (Pfa design).
%
%   [nInt, nRef] = rsp_cfar_looks(P, replica)
%
%   P       : chain parameters
%   replica : matched-filter taps that shape the range noise (long pulse),
%             [] = matched filter off (white noise in range)
%   nInt    : 1 x nfft, effective looks of one integrated cell
%   nRef    : 1 x nfft, effective independent cells of one reference window
%
%   Every block is linear up to the detector, so the correlation of the
%   noise samples summed by the integrator and by the CFAR window is known:
%     range  : autocorrelation of the matched filter
%     frames : one FFT frame of bin b is sum_p g_b(p) x(p + start) with
%              g_b = canceler taps (*) (window .* bin phasor); frames are
%              hop pulses apart
%   A sum of squared correlated Gaussian samples is sum_i lambda_i E_i
%   (lambda = eigenvalues of their correlation matrix, E_i exponential).
%     nInt : Gamma shape with the same mean and variance as the integrated cell
%     nRef : Gamma shape whose Laplace transform equals the exact one of the
%            reference window at the operating threshold. This is what sets
%            the false-alarm rate (for CA with nInt = 1 it is exact), and it
%            is solved jointly with the threshold factor of P.cfar.pfa.
%   P.cfar.nIntEffective / P.cfar.nRefEffective override these values.

C    = P.cfar;
N    = C.nRef;
n    = P.nci.nFrames;
nfft = max(P.fft.nfft, P.fft.nPulses);

% --- range correlation of the reference cells -----------------------------------
if isempty(replica) || ~P.mf.enable
    lamR = ones(N, 1);
else
    h    = replica(:);
    acf  = conv(h, conj(flipud(h)));
    acf  = acf(numel(h):end) / acf(numel(h));          % rho at lags 0, 1, 2 ...
    rho  = zeros(N, 1);
    m    = min(N, numel(acf));
    rho(1:m) = acf(1:m);
    lamR = hermEig(rho);
end

% --- frame correlation of the integrated cell, per bin -----------------------------
if P.canceler.enable
    c = rsp_canceler_taps(P.canceler);
else
    c = 1;
end
w    = rsp_window(P.fft.window, P.fft.nPulses, P.fft.windowParam);
p    = (0:P.fft.nPulses-1)';
k    = (0:nfft-1) - P.fft.shift * floor(nfft/2);     % frequency index of each bin
lamF = zeros(n, nfft);
for b = 1:nfft
    g   = conv(w .* exp(-1j*2*pi*k(b)*p/nfft), c(:));
    acf = conv(g, conj(flipud(g)));
    acf = acf(numel(g):end) / acf(numel(g));
    lag = (0:n-1)' * P.fft.hop;
    rho = zeros(n, 1);
    in  = lag < numel(acf);
    rho(in) = acf(lag(in) + 1);
    lamF(:, b) = hermEig(rho);
end

% --- effective looks ------------------------------------------------------------------
nInt = sum(lamF, 1).^2 ./ sum(lamF.^2, 1);
if ~isempty(C.nIntEffective)
    nInt(:) = C.nIntEffective;
end
nRef = zeros(1, nfft);
[~, first, same] = unique(round(lamF' * 1e8) / 1e8, 'rows');   % bins with equal statistics
for u = 1:numel(first)
    b = first(u);
    e = kron(lamR, lamF(:, b));                          % eigenvalues of one window
    if ~isempty(C.nRefEffective)
        v = C.nRefEffective;
    elseif strcmpi(C.thresholdMode, 'pfa') && strcmpi(P.nci.law, 'square')
        v = laplaceMatch(e, nInt(b), C);
    else
        v = sum(e)^2 / sum(e.^2) / nInt(b);              % variance match
    end
    nRef(same == u) = v;
end
end

function nRef = laplaceMatch(e, K, C)
% Shape a of a Gamma window sum with mean m = sum(e) whose Laplace transform
% equals the exact one, prod(1 + s e)^-1, at the threshold point
% s = alpha*K/m (decay rate of the false-alarm integrand). alpha depends on
% a, so both are iterated to a fixed point.
m    = sum(e);
a    = m^2 / sum(e.^2);                                  % start: variance match
for it = 1:20
    alpha = rsp_cfar_factor(C.type, a / K, K, C.pfa);
    s     = alpha * K / m;
    T     = sum(log1p(s * e));                            % -log of exact transform
    f     = @(x) exp(x) .* log1p(s * m ./ exp(x)) - T;    % gamma(-log transform) - T
    aNew  = exp(fzero(f, [log(1e-3) log(1e3 * m)]));
    if abs(aNew - a) < 1e-4 * a
        a = aNew;
        break
    end
    a = aNew;
end
nRef = a / K;
end

function lam = hermEig(rho)
% Eigenvalues of the Hermitian Toeplitz correlation matrix with first column rho.
rho(1) = real(rho(1));
R   = toeplitz(rho, conj(rho));
lam = max(real(eig((R + R') / 2)), 0);
end
