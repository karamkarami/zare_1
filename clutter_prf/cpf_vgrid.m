function v = cpf_vgrid(cfg)
%CPF_VGRID Target radial velocities of the criteria: +-[vMin vMax] [m/s] (row).

p = cfg.target.vMinMps:cfg.prf.dvMps:cfg.target.vMaxMps;
if p(end) < cfg.target.vMaxMps
    p(end+1) = cfg.target.vMaxMps;
end
v = [-fliplr(p) p];
end
