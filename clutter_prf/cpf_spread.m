function [sScanHz, kMot] = cpf_spread(st)
%CPF_SPREAD Spectrum spread of the ground clutter caused by the radar itself.
%
%   [sScanHz, kMot] = cpf_spread(st)
%
%   sScanHz : sigma of the antenna scan modulation [Hz] = 0.265 * omega / beamwidth
%             (Gaussian two-way beam; omega measured from the azimuths, else rpm)
%   kMot    : platform motion spread per m/s of velocity across the beam [Hz/(m/s)]
%             sigma_mot = kMot * |v x u| = 2/lambda * 0.300 * beamwidth[rad] * v_perp
%             (Doppler 2 v cos(az - heading)/lambda changes across the beam)

cfg = st.cfg;
omega = st.omegaDps;
if isempty(omega)
    omega = 6 * cfg.radar.rpm;
end
sScanHz = 0.265 * omega / cfg.radar.beamwidthDeg;
kMot = 2 / st.lambda * 0.300 * cfg.radar.beamwidthDeg * pi / 180;
end
