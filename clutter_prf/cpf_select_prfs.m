function sel = cpf_select_prfs(cfg, design, lambda)
%CPF_SELECT_PRFS Choose the P PRFs of the scan-to-scan set.
%
%   sel = cpf_select_prfs(cfg, design, lambda)
%
%   Candidates: cfg.prf.minHz : stepHz : maxHz. For every candidate the filter
%   is designed for the clutter in design (.cnrDb, .sigmaHz), then the set of
%   cfg.prf.count PRFs with the lowest cpf_set_cost is searched:
%     - worst-case clutter attenuation (worst PRF of the set) against the
%       attenuation needed (cpf_ca_required)
%     - ripple of the mean response over +-[vMin vMax]
%     - at most cfg.prf.maxBlind PRFs blind at any velocity (cfg.prf.blindMode)
%     - cumulative Pd over the P scans (reported, in the cost with wPd > 0)
%   Every subset is tried when there are at most cfg.prf.maxExhaustive of them,
%   otherwise a local search (exchange one PRF at a time) from several starts.
%
%   sel: .prfHz (ascending), .caDb, .rippleDb, .caWorstDb, .caReqDb, .blindMax,
%        .pcMin, .pcMean, .J, .feasible, .v, .G (P x nV gains), .W (taps), .design

cand = cfg.prf.minHz:cfg.prf.stepHz:cfg.prf.maxHz;
K = numel(cand);
P = cfg.prf.count;
if P > K
    error('cpf_select_prfs:count', 'Only %d candidate PRFs for P = %d.', K, P);
end
M = cpf_prf_metrics(cfg, design, cand, lambda);
nV = numel(M.v);
caReq = cpf_ca_required(cfg, design);
cost = @(S) cpf_set_cost(cfg, permute(reshape(M.G(S.', :), [P size(S, 1) nV]), [2 1 3]), ...
                       permute(reshape(M.blind(S.', :), [P size(S, 1) nV]), [2 1 3]), ...
                       reshape(M.caDb(S), size(S)), caReq);

nSub = exp(gammaln(K + 1) - gammaln(P + 1) - gammaln(K - P + 1));
if nSub <= cfg.prf.maxExhaustive
    all = nchoosek(1:K, P);
    best = [];  bestJ = Inf;
    for i0 = 1:500:size(all, 1)
        S = all(i0:min(i0 + 499, end), :);
        c = costOf(cost, S);
        [j, i] = min(c);
        if j < bestJ
            bestJ = j;  best = S(i, :);
        end
    end
else
    seed = 12345;
    best = [];  bestJ = Inf;
    for r = 1:max(cfg.prf.restarts, 1)
        if r == 1
            cur = unique(round(linspace(1, K, P)));
        else
            [p, seed] = lcgPerm(K, seed);
            cur = sort(p(1:P));
        end
        if numel(cur) < P
            continue
        end
        curJ = costOf(cost, cur);
        improved = true;
        while improved
            improved = false;
            for i = 1:P
                others = setdiff(1:K, cur);
                S = repmat(cur, numel(others), 1);
                S(:, i) = others(:);
                c = costOf(cost, S);
                [j, k] = min(c);
                if j < curJ - 1e-9
                    cur = sort(S(k, :));  curJ = j;  improved = true;
                end
            end
        end
        if curJ < bestJ
            bestJ = curJ;  best = cur;
        end
    end
end
sel = cpf_prf_eval(cfg, design, cand(sort(best)), lambda);
sel.candidates = cand;
sel.searched   = min(nSub, cfg.prf.maxExhaustive);
end

function J = costOf(cost, S)
% S: one PRF set (candidate indices) per row
c = cost(S);
J = c.J;
end

function [p, seed] = lcgPerm(K, seed)
% Reproducible random permutation (does not touch the global random state).
u = zeros(1, K);
for i = 1:K
    seed = mod(1103515245 * seed + 12345, 2^31);
    u(i) = seed / 2^31;
end
[~, p] = sort(u);
end
