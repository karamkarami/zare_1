function m = cpf_map_empty(cfg, nA, nR)
%CPF_MAP_EMPTY Clutter map with no data: nA azimuth x nR range cells.
%
%   Per cell (index ia + (ir - 1) * nA): ground-frame clutter of the cell
%     C    clutter power per sample (linear, same unit as the data)
%     sig2 intrinsic spectrum spread sigma_v^2 [m^2/s^2] (no scan or platform spread)
%     vm   mean velocity [m/s] (0 for fixed ground)
%     n    number of updates, nw of them with a width measurement

m.nA = nA;  m.nR = nR;
m.C    = zeros(nA * nR, 1);
m.sig2 = cfg.clutter.sigmaVMps^2 * ones(nA * nR, 1);
m.vm   = zeros(nA * nR, 1);
m.n    = zeros(nA * nR, 1);
m.nw   = zeros(nA * nR, 1);
m.origin = [];                          % ground point of the map centre [east north]
end
