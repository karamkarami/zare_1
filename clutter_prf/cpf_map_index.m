function L = cpf_map_index(st, pos, azDeg, rangeM)
%CPF_MAP_INDEX Map cell of every sample (pulses x range), 0 outside the map.
%
%   The ground point of a sample is the platform position plus the range along
%   the antenna azimuth; its map cell is the polar cell of that point seen from
%   the map origin. For a stationary radar at the origin this is simply the
%   (azimuth, range) cell.

cfg = st.cfg;
o  = st.map.origin;
E  = (pos(:, 1) - o(1)) + sind(azDeg(:)) * rangeM(:).';
N  = (pos(:, 2) - o(2)) + cosd(azDeg(:)) * rangeM(:).';
a  = mod(atan2(E, N) * 180 / pi, 360);
ia = min(floor(a / cfg.map.azCellDeg) + 1, st.map.nA);
ir = floor(hypot(E, N) / cfg.map.rangeCellM) + 1;
L  = ia + (ir - 1) * st.map.nA;
L(ir > st.map.nR) = 0;
end
