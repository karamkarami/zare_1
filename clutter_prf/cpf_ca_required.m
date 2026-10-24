function ca = cpf_ca_required(cfg, design)
%CPF_CA_REQUIRED Clutter attenuation the filters need [dB].
%
%   cfg.prf.caMinDb, or the design CNR when it is empty: the clutter is then
%   brought down to the noise.

if isempty(cfg.prf.caMinDb)
    ca = design.cnrDb;
else
    ca = cfg.prf.caMinDb;
end
end
