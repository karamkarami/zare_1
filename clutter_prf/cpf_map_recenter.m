function st = cpf_map_recenter(st, origin)
%CPF_MAP_RECENTER Move the map origin to a new ground point (moving platform).
%
%   Every cell of the new grid takes the old cell under its centre, so the
%   clutter learned so far is kept.

cfg = st.cfg;
m = st.map;
[ia, ir] = ndgrid(1:m.nA, 1:m.nR);
a = (ia(:) - 0.5) * cfg.map.azCellDeg;
r = (ir(:) - 0.5) * cfg.map.rangeCellM;
E = origin(1) - m.origin(1) + r .* sind(a);
N = origin(2) - m.origin(2) + r .* cosd(a);
ja = min(floor(mod(atan2(E, N) * 180 / pi, 360) / cfg.map.azCellDeg) + 1, m.nA);
jr = floor(hypot(E, N) / cfg.map.rangeCellM) + 1;
old = ja + (jr - 1) * m.nA;
in  = jr <= m.nR;
n = cpf_map_empty(cfg, m.nA, m.nR);
for f = {'C', 'sig2', 'vm', 'n', 'nw'}
    n.(f{1})(in) = m.(f{1})(old(in));
end
n.origin = origin;
st.map = n;
end
