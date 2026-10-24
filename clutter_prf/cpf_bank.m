function bank = cpf_bank(cfg, prfHz, lambda, need)
%CPF_BANK Filters of one PRF over the grid of clutter width and CNR.
%
%   bank = cpf_bank(cfg, prfHz, lambda)          empty bank of this PRF
%   bank = cpf_bank(cfg, bank, lambda, need)     design the entries in need
%                                                that are not designed yet
%
%   bank.W(:, b)  taps for width cfg.mti.sGrid(iS) and CNR cfg.mti.cnrGridDb(iC),
%                 b = iS + (iC - 1) * numel(sGrid); designed on first use
%   bank.caDb(b)  clutter attenuation of that filter
%   b = 0 in a cell means no filter (cpf_process).

sG = cfg.mti.sGrid;
cG = cfg.mti.cnrGridDb;
nS = numel(sG);
if nargin < 4
    bank.prfHz = prfHz;
    bank.B     = cpf_band(cfg, prfHz, lambda);
    bank.W     = zeros(cfg.mti.maxLength, nS * numel(cG));
    bank.caDb  = nan(1, nS * numel(cG));
    bank.done  = false(1, nS * numel(cG));
    return
end
bank = prfHz;
need = unique(need(need > 0));
for b = need(~bank.done(need))
    iS = mod(b - 1, nS) + 1;
    iC = (b - iS) / nS + 1;
    [bank.W(:, b), bank.caDb(b)] = cpf_design_filter(bank.B, sG(iS), cG(iC), cfg);
    bank.done(b) = true;
end
end
