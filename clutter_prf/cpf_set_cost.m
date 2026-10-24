function c = cpf_set_cost(cfg, G, blind, caDb, caReqDb)
%CPF_SET_COST Criteria and cost of PRF sets.
%
%   c = cpf_set_cost(cfg, G, caDb)
%
%   G    : nSets x P x nV filter gains (linear) of the PRFs of each set
%   blind: nSets x P x nV, true where the target is blind (cpf_prf_metrics)
%   caDb : nSets x P clutter attenuation of those PRFs
%   caReqDb : clutter attenuation needed (clutter brought down to the noise)
%
%   c.rippleDb : ripple of the mean response (mean over the P PRFs) over the band
%   c.caWorstDb: clutter attenuation of the worst PRF
%   c.blindMax : largest number of PRFs in the clutter notch at one velocity
%   c.overPct  : % of the band where more than cfg.prf.maxBlind PRFs are blind
%   c.pcMin    : smallest cumulative Pd over the P scans (Swerling 1 scan to scan,
%                SNR scaled by the filter gain, 0 where blind)
%   c.pcMean   : mean cumulative Pd over the band
%   c.J        : cost (lower is better):
%                wRipple * rippleDb + wCa * (dB the worst PRF misses caReqDb by)
%                - wPd * 100 * pcMin
%                + 1000 per PRF over maxBlind + 10 per % of band over maxBlind
%                (when maxBlind cannot be met, the set with the smallest
%                over-limit part of the band wins)
%                CA beyond caReqDb is not rewarded: the clutter is then already
%                below the noise.

[nSets, P, nV] = size(G);
blindIn = blind;
Gm = reshape(mean(G, 2), nSets, nV);
c.rippleDb  = 10 * log10(max(Gm, [], 2) ./ max(min(Gm, [], 2), 1e-30));
c.caWorstDb = min(caDb, [], 2);
blind       = reshape(sum(blind, 2), nSets, nV);
c.blindMax  = max(blind, [], 2);
c.overPct   = 100 * mean(blind > cfg.prf.maxBlind, 2);
snr  = 10^(cfg.prf.snrDb / 10);
pd   = cfg.prf.pfa .^ (1 ./ (1 + snr * G));
pd(blindIn) = 0;
pc   = 1 - reshape(prod(1 - pd, 2), nSets, nV);
c.pcMin  = min(pc, [], 2);
c.pcMean = mean(pc, 2);
c.caShortDb = max(caReqDb - c.caWorstDb, 0);
c.J = cfg.prf.wRipple * c.rippleDb + cfg.prf.wCa * c.caShortDb - cfg.prf.wPd * 100 * c.pcMin ...
    + 1000 * max(c.blindMax - cfg.prf.maxBlind, 0) + 10 * c.overPct;
end
