function [prf, st, sel] = cpf_next_prf(st, speedMps)
%CPF_NEXT_PRF PRF of the next scan; re-selects the PRF set from time to time.
%
%   [prf, st, sel] = cpf_next_prf(st)            speed = last platform speed in the data
%   [prf, st, sel] = cpf_next_prf(st, speedMps)
%
%   Call once before every scan. The PRFs of st.prfSet are used in turn.
%   Every cfg.update.everyScans scans, or when the platform speed changed by
%   cfg.update.speedDeltaMps since the last selection, the set is selected
%   again for the clutter now in the map (cpf_design_clutter). The new set
%   replaces the old one only if its cost is lower by cfg.update.hysteresis,
%   or if the old set no longer meets its limits. The map is kept either way:
%   it is stored in m/s, so the new PRFs are filtered right away.
%   sel: the new selection when the set was re-evaluated, else [].

cfg = st.cfg;
if nargin < 2 || isempty(speedMps)
    speedMps = norm(st.vel);
end
st.schedScan = st.schedScan + 1;
sel = [];
due = st.schedScan - st.selScan > cfg.update.everyScans || ...
      abs(speedMps - st.selSpeed) >= cfg.update.speedDeltaMps;
if due
    d   = cpf_design_clutter(st, speedMps);
    new = cpf_select_prfs(cfg, d, st.lambda);
    cur = cpf_prf_eval(cfg, d, st.prfSet, st.lambda);
    if new.J < cur.J - cfg.update.hysteresis || ~cur.feasible && new.J < cur.J
        st.prfSet  = new.prfHz;
        st.sel     = new;
        st.prfNext = 1;
        new.changed = true;
    else
        st.sel = cur;
        new = cur;
        new.changed = false;
    end
    sel = new;
    st.selScan  = st.schedScan;
    st.selSpeed = speedMps;
end
prf = st.prfSet(st.prfNext);
st.prfNext = mod(st.prfNext, numel(st.prfSet)) + 1;
end
