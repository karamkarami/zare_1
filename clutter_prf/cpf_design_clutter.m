function d = cpf_design_clutter(st, speedMps)
%CPF_DESIGN_CLUTTER Clutter the PRF set is designed for.
%
%   d = cpf_design_clutter(st, speedMps)
%
%   d.cnrDb   : cfg.update.cnrPercentile of the map CNR (cells with clutter),
%               cfg.clutter.cnrDb while the map is empty; plus cfg.mti.cnrMarginDb
%   d.sigmaHz : total spread for the PRF design [Hz]:
%               widthMargin * sqrt(intrinsic^2 + scan^2 + platform^2), the
%               intrinsic part from the map (cfg.update.sigmaPercentile), the
%               platform part at speedMps, broadside ('worst') or averaged over
%               azimuth ('mean')

cfg = st.cfg;
m = st.map;
cnrDb = cfg.clutter.cnrDb;
sig2  = cfg.clutter.sigmaVMps^2;
if ~isempty(st.noise) && st.noise > 0
    ok = m.n >= cfg.map.minUpdates & m.C > st.noise * 10^(cfg.mti.bypassCnrDb / 10);
    if nnz(ok) >= 10
        cnrDb = pct(10 * log10(m.C(ok) / st.noise), cfg.update.cnrPercentile);
        okw = ok & m.nw > 0;
        if nnz(okw) >= 10
            sig2 = pct(m.sig2(okw), cfg.update.sigmaPercentile);
        end
    end
end
[sScan, kMot] = cpf_spread(st);
if strcmpi(cfg.update.motionCase, 'mean')
    vPerp2 = speedMps^2 / 2;            % mean of sin^2 over azimuth
else
    vPerp2 = speedMps^2;                % broadside
end
d.cnrDb   = cnrDb + cfg.mti.cnrMarginDb;
d.sigmaHz = cfg.mti.widthMargin * sqrt((2 / st.lambda)^2 * sig2 + sScan^2 + kMot^2 * vPerp2);
d.speedMps = speedMps;
end

function p = pct(x, q)
% percentile without a toolbox (nearest rank)
x = sort(x(:));
p = x(min(max(ceil(q / 100 * numel(x)), 1), numel(x)));
end
