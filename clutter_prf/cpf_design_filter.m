function [w, caDb, edge] = cpf_design_filter(B, s, cnrDb, cfg)
%CPF_DESIGN_FILTER Slow-time FIR for one PRF, clutter width and CNR.
%
%   [w, caDb, edge] = cpf_design_filter(B, s, cnrDb, cfg)
%
%   B      : target band of the PRF (cpf_band)
%   s      : clutter spectrum sigma / PRF (Gaussian, centred on 0 Hz)
%   cnrDb  : clutter-to-noise ratio per sample [dB]
%
%   w      : FIR taps (y(m) = sum_k w(k+1) x(m-k)), maxLength x 1, unit norm
%            (unit noise gain); an 'auto' length filter is padded with zeros
%   caDb   : clutter attenuation 1 / (w' Rc w), Rc = unit-power clutter correlation
%   edge   : clutter edge / PRF, where the clutter spectrum falls to the noise
%
%   1. Clutter: the taps are kept in the span of the eigenvectors of Rc whose
%      clutter power CNR * eig stays under cfg.mti.residueDb (eigencanceler),
%      so any unit-norm filter there leaves at most that residue.
%   2. Passband (target velocities whose folded Doppler is outside the clutter
%      edge): inside that span the taps maximise the mean gain in dB (the
%      geometric mean of the gain). A low gain anywhere costs a lot, so the
%      notch stays narrow and the passband flat. Maximising the plain mean gain
%      would give the widest notch instead, because for unit-norm taps the mean
%      gain over all Doppler is always 1.

Nmax = size(B.S, 2);
cnr  = 10^(cnrDb / 10);
edge = s * sqrt(2 * log(max(cnr / (sqrt(2 * pi) * s), 1)));
pass = B.dn >= edge;
if nnz(pass) < Nmax
    pass = true(size(B.dn));          % clutter wider than the band: use it all
end
rc = exp(-2 * pi^2 * s^2 * (0:Nmax-1).^2);
allow = 10^(cfg.mti.residueDb / 10);
if strcmpi(cfg.mti.lengthMode, 'auto')
    lengths = 2:Nmax;
else
    lengths = Nmax;
end
for N = lengths
    Rc = toeplitz(rc(1:N));
    [E, lam] = eig((Rc + Rc') / 2);
    [lam, i] = sort(real(diag(lam)));
    E = E(:, i);
    D = max(nnz(cnr * lam <= allow), 1);
    if cnr * lam(1) <= allow && D >= min(2, N - 1)
        break                         % this length meets the residue with room to shape
    end
end
U = E(:, 1:D);
C = B.S(pass, 1:N) * U;               % gain at velocity v: |C(v,:) * z|^2, w = U z
z = meanDbGain(C, cfg.mti.designIter);
w = zeros(Nmax, 1);
w(1:N) = U * z;
w = w / norm(w);
if max(abs(imag(w))) < 1e-9 * max(abs(w))
    w = real(w);
end
caDb = -10 * log10(max(real(w(1:N)' * Rc * w(1:N)), 1e-30));
end

function z = meanDbGain(C, nIter)
% Maximise L(z) = mean(log g), g = |C z|^2, |z| = 1, by projected gradient
% ascent with step control, starting from the best mean (linear) gain.
D = size(C, 2);
[V, d] = eig(C' * C);
[~, i] = max(real(diag(d)));
z = V(:, i);
if D == 1
    return
end
[L, grad] = objective(C, z);
mu = 0.5;
for it = 1:nIter
    d = grad - (z' * grad) * z;                % tangent to the unit sphere
    nd = norm(d);
    if nd < 1e-10
        break
    end
    while mu > 1e-8
        zn = z + mu * d / nd;
        zn = zn / norm(zn);
        [Ln, gn] = objective(C, zn);
        if Ln > L
            break
        end
        mu = mu / 2;
    end
    if mu <= 1e-8
        break
    end
    gain = Ln - L;
    z = zn;  L = Ln;  grad = gn;
    mu = min(mu * 1.5, 1);
    if gain < 1e-7
        break
    end
end
end

function [L, grad] = objective(C, z)
y = C * z;
g = max(real(y).^2 + imag(y).^2, 1e-15);
n = numel(g);
L = sum(log(g)) / n;
grad = C' * (y ./ g) / n;
end
