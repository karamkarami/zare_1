function B = cpf_band(cfg, prfHz, lambda)
%CPF_BAND Target band of one PRF: steering vectors and folded distance to the clutter.
%
%   B = cpf_band(cfg, prfHz, lambda)
%
%   B.v  : velocities +-[vMin vMax] over ground [m/s] (after motion compensation
%          the ground clutter is at 0 Hz, so a target is at 2 v / lambda from it)
%   B.S  : steering vectors, numel(v) x maxLength (cpf_steering)
%   B.dn : distance of the folded Doppler to the nearest clutter line k*PRF,
%          divided by the PRF (0 ... 0.5)

B.prfHz = prfHz;
B.v  = cpf_vgrid(cfg);
f    = 2 * B.v(:) / lambda;
B.S  = cpf_steering(f, prfHz, cfg.mti.maxLength);
B.dn = abs(mod(f / prfHz + 0.5, 1) - 0.5);
end
