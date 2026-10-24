function [x, y] = rsp_clutter_outline(L)
%RSP_CLUTTER_OUTLINE Outline of a clutter layer for drawing, x = east, y = north [km].
%
%   [x, y] = rsp_clutter_outline(L)      L : one layer (rsp_clutter_types)

w = L.azDeg(2) - L.azDeg(1);
if abs(w) >= 360
    a = linspace(0, 360, 181);
else
    a = L.azDeg(1) + linspace(0, mod(w, 360), max(2, ceil(mod(w, 360)) + 1));
end
r1 = L.rangeM(1) / 1e3;
r2 = L.rangeM(2) / 1e3;
if abs(w) >= 360
    x = [r2 * sind(a), NaN, r1 * sind(a)];
    y = [r2 * cosd(a), NaN, r1 * cosd(a)];
else
    x = [r1 * sind(a), r2 * sind(fliplr(a)), r1 * sind(a(1))];
    y = [r1 * cosd(a), r2 * cosd(fliplr(a)), r1 * cosd(a(1))];
end
end
