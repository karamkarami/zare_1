function [alpha, pfaFun] = rsp_cfar_factor(type, nRef, nInt, pfa)
%RSP_CFAR_FACTOR Threshold multiplier of a CA / SO / GO CFAR for a target Pfa.
%
%   [alpha, pfaFun] = rsp_cfar_factor(type, nRef, nInt, pfa)
%
%   type   : 'SO' | 'CA' | 'GO'
%   nRef   : (effective) reference cells on each side, may be non-integer
%   nInt   : (effective) looks integrated non-coherently, may be non-integer
%   pfa    : desired probability of false alarm
%   alpha  : multiplier applied to the noise estimate (mean of the cells)
%   pfaFun : @(alpha) Pfa, for checking other multipliers
%
%   Model: square-law detector, Gaussian noise. The cell under test is
%   Gamma(nInt) and one reference window sum is Gamma(a), a = nRef*nInt.
%   With Z the selected window sum (one window for SO/GO, both for CA):
%       Pfa(alpha) = E_Z[ Q(nInt, alpha*Z/nRef) ]
%   with Q the regularised upper incomplete gamma function. The expectation
%   is a fixed quadrature on a logarithmic grid (fast and accurate for very
%   large multipliers); alpha is solved with fzero.
%   For CA and nInt = 1 this is exactly (1 + alpha/(2 nRef))^(-2 nRef), for
%   SO it reproduces the closed form of Weiss (1982).

type = upper(type);
K    = nInt;
a    = nRef * K;                                  % shape of one window sum

switch type
    case {'SO', 'GO'}
        sc = nRef;
        sh = a;
    case 'CA'
        sc = 2*nRef;
        sh = 2*a;
    otherwise
        error('rsp_cfar_factor:type', 'Unknown CFAR type "%s".', type);
end

% log-spaced grid over the support of Z; integrand f(z) dz = f(z) z dlog(z)
u  = linspace(log(sh) - 40, log(sh + 40*sqrt(sh) + 50), 6000)';
z  = exp(u);
du = u(2) - u(1);
switch type
    case 'SO'
        F  = gammainc(z, a);
        fz = 2*exp(a*log(z) - z - gammaln(a)).*(1 - F);
        p0 = 1 - (1 - F(1))^2;                    % mass below the grid (Q ~ 1)
    case 'GO'
        F  = gammainc(z, a);
        fz = 2*exp(a*log(z) - z - gammaln(a)).*F;
        p0 = F(1)^2;
    case 'CA'
        fz = exp(sh*log(z) - z - gammaln(sh));
        p0 = gammainc(z(1), sh);
end
wq = fz * du;
wq([1 end]) = wq([1 end]) / 2;                    % trapezoid weights

pfaFun = @(al) p0 + gammainc(al*z/sc, K, 'upper')' * wq;

if nargin < 4 || isempty(pfa)
    alpha = [];
    return
end
h     = @(x) log(max(pfaFun(exp(x)), realmin)) - log(pfa);   % finite at both ends
alpha = exp(fzero(h, [log(1e-3) log(1e12)]));
end
