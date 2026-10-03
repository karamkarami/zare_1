function out = rsp_cfar(S, C, law, nInt, nRefEff)
%RSP_CFAR Range CFAR (smallest-of by default) on every frame and Doppler bin.
%
%   out = rsp_cfar(S, C, law, nInt, nRefEff)
%
%   S    : frames x range x doppler, integrated detector output (rsp_nci)
%   C    : P.cfar (see rsp_default_params)
%   law  : detector law used by the integrator ('square' | 'linear' | 'log')
%   nInt    : effective looks in S, scalar or one per Doppler bin (default 1)
%   nRefEff : effective independent cells per reference window, scalar or
%             one per Doppler bin (default C.nRef)
%             nInt and nRefEff are only used by thresholdMode = 'pfa',
%             see rsp_cfar_looks
%
%   out.det        logical, detections                      (size of S)
%   out.threshold  single, detection threshold              (size of S)
%   out.map        single, S on detections and 0 elsewhere  (size of S)
%   out.factor     threshold multiplier of each bin (dB offset for the 'log' law)
%   out.list       table-like struct of detections: frame, range, bin,
%                  value, threshold, snrDb (value over noise estimate)
%
%   For cell under test r the two reference windows are
%       lead : r-G-N .. r-G-1        lag : r+G+1 .. r+G+N
%   (N = C.nRef, G = C.nGuard). Their means are computed for all cells at
%   once with a cumulative sum, so the cost does not depend on N.
%       SO : noise = min(lead, lag)    CA : (lead + lag)/2    GO : max(lead, lag)
%   Threshold = factor * noise  (noise + factor for the 'log' law).

if nargin < 3 || isempty(law),     law     = 'square'; end
if nargin < 4 || isempty(nInt),    nInt    = 1;        end
if nargin < 5 || isempty(nRefEff), nRefEff = C.nRef;   end

[nF, R, nB] = size(S);
N    = C.nRef;
G    = C.nGuard;
type = upper(C.type);
isLog = strcmpi(law, 'log');

% --- threshold multiplier ----------------------------------------------------
switch lower(C.thresholdMode)
    case 'pfa'
        if ~strcmpi(law, 'square')
            error('rsp_cfar:pfa', ['thresholdMode = ''pfa'' needs P.nci.law = ''square''. ' ...
                'Use thresholdMode = ''factor'' for the ''%s'' law.'], law);
        end
        looks = round([nInt(:) .* ones(nB, 1), nRefEff(:) .* ones(nB, 1)] * 1e4) / 1e4;
        [u, ~, j] = unique(looks, 'rows');              % one solve per distinct pair
        factor = zeros(1, nB);
        for i = 1:size(u, 1)
            factor(j == i) = rsp_cfar_factor(type, u(i, 2), u(i, 1), C.pfa);
        end
    case 'factor'
        if isLog
            factor = C.factorDb;                         % additive, in dB
        elseif strcmpi(law, 'linear')
            factor = 10^(C.factorDb/20);
        else
            factor = 10^(C.factorDb/10);
        end
        factor = repmat(factor, 1, nB);
    otherwise
        error('rsp_cfar:mode', 'Unknown P.cfar.thresholdMode "%s".', C.thresholdMode);
end

% --- reference window bounds (computed once, shared by all rows) ------------
r      = 1:R;
leadLo = max(r - G - N, 1);   leadHi = min(r - G - 1, R);
lagLo  = max(r + G + 1, 1);   lagHi  = min(r + G + N, R);
nLead  = max(leadHi - leadLo + 1, 0);
nLag   = max(lagHi  - lagLo  + 1, 0);
leadHi = max(leadHi, leadLo - 1);              % empty window -> zero sum
lagLo  = min(lagLo,  lagHi + 1);

fullLead = nLead == N;
fullLag  = nLag  == N;
switch lower(C.edge)
    case 'onesided'
        useLead = fullLead | (~fullLag & nLead > 0);
        useLag  = fullLag  | (~fullLead & nLag > 0);
    case 'partial'
        useLead = nLead > 0;
        useLag  = nLag  > 0;
    case 'none'
        useLead = fullLead & fullLag;
        useLag  = useLead;
    otherwise
        error('rsp_cfar:edge', 'Unknown P.cfar.edge "%s".', C.edge);
end
tested = useLead | useLag;
if ~isempty(C.rangeCells)
    tested = tested & r >= C.rangeCells(1) & r <= C.rangeCells(end);
end
bins = C.bins;
if isempty(bins)
    bins = 1:nB;
end

% --- CFAR, one Doppler bin at a time (bounded memory) ------------------------
out.det       = false(nF, R, nB);
out.threshold = inf(nF, R, nB, 'single');            % Inf = cell not tested
out.factor    = factor;
for b = bins
    x  = double(S(:, :, b));
    cs = [zeros(nF, 1) cumsum(x, 2)];
    mLead = (cs(:, leadHi + 1) - cs(:, leadLo)) ./ nLead;
    mLag  = (cs(:, lagHi  + 1) - cs(:, lagLo))  ./ nLag;
    mLead(:, ~useLead) = NaN;
    mLag(:,  ~useLag)  = NaN;
    switch type
        case 'SO', noise = min(mLead, mLag);          % NaN-aware: one side only
        case 'GO', noise = max(mLead, mLag);
        case 'CA', noise = meanNan(mLead, mLag);
        otherwise, error('rsp_cfar:type', 'Unknown P.cfar.type "%s".', C.type);
    end
    if isLog
        thr = noise + factor(b);
    else
        thr = noise * factor(b);
    end
    thr(:, ~tested) = Inf;
    out.det(:, :, b)       = x > thr;
    out.threshold(:, :, b) = thr;
end

out.map = zeros(nF, R, nB, 'single');
out.map(out.det) = S(out.det);

% --- detection list -----------------------------------------------------------
k = find(out.det);
[f, rr, bb] = ind2sub([nF R nB], k);
val   = double(S(k));
thr   = double(out.threshold(k));
fb    = factor(bb);
fb    = fb(:);
if isLog
    snr = val - (thr - fb);
else
    snr = 10*log10(val ./ (thr ./ fb)) * (1 + strcmpi(law, 'linear'));
end
out.list = struct('frame', f, 'range', rr, 'bin', bb, 'value', val, ...
                  'threshold', thr, 'snrDb', snr);
end

function m = meanNan(a, b)
% Mean of two arrays, ignoring NaN entries.
na = ~isnan(a);
nb = ~isnan(b);
a(~na) = 0;
b(~nb) = 0;
m = (a + b) ./ (na + nb);
m(~(na | nb)) = NaN;
end
