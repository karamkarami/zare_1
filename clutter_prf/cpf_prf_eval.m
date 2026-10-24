function sel = cpf_prf_eval(cfg, design, prfHz, lambda)
%CPF_PRF_EVAL Criteria of a given PRF set (same fields as cpf_select_prfs).
%
%   sel = cpf_prf_eval(cfg, design, prfHz, lambda)

M = cpf_prf_metrics(cfg, design, prfHz, lambda);
P = numel(prfHz);
c = cpf_set_cost(cfg, reshape(M.G, [1 P numel(M.v)]), reshape(M.blind, [1 P numel(M.v)]), ...
                 M.caDb, cpf_ca_required(cfg, design));
sel = c;
sel.prfHz  = M.prfHz;
sel.caDb   = M.caDb;
sel.v      = M.v;
sel.G      = M.G;
sel.blind  = M.blind;
sel.W      = M.W;
sel.design = design;
sel.design.fcHz = 299792458 / lambda;
sel.maxBlind = cfg.prf.maxBlind;
sel.caReqDb  = cpf_ca_required(cfg, design);
Gm = mean(M.G, 1);
[~, i] = min(Gm);
sel.vWorstMean = M.v(i);                 % velocity of the lowest mean response
pd = cfg.prf.pfa .^ (1 ./ (1 + 10^(cfg.prf.snrDb / 10) * M.G));
pd(M.blind) = 0;
pc = 1 - prod(1 - pd, 1);
[~, i] = min(pc);
sel.vWorstPd = M.v(i);                   % velocity of the lowest cumulative Pd
sel.feasible = c.blindMax <= cfg.prf.maxBlind && c.caShortDb == 0;
end
