function M = cpf_prf_metrics(cfg, design, prfHz, lambda)
%CPF_PRF_METRICS Filter, clutter attenuation and velocity response of each PRF.
%
%   M = cpf_prf_metrics(cfg, design, prfHz, lambda)
%
%   design : .cnrDb and .sigmaHz of the clutter to design for (cpf_design_clutter)
%   M.v    : velocity grid +-[vMin vMax] [m/s]
%   M.G    : numel(prfHz) x numel(v) filter power gain (linear, re noise gain)
%   M.blind: numel(prfHz) x numel(v), true where the target is blind at that PRF:
%            cfg.prf.blindMode 'clutter': its folded Doppler is inside the
%            clutter region (clutter edge + half a Doppler FFT bin);
%            'notch': the filter gain is below cfg.prf.blindDb
%   M.caDb : clutter attenuation of each PRF's filter [dB]
%   M.W    : taps, maxLength x numel(prfHz)

M.v     = cpf_vgrid(cfg);
M.prfHz = prfHz(:).';
K = numel(prfHz);
M.G    = zeros(K, numel(M.v));
M.blind = false(K, numel(M.v));
M.caDb = zeros(1, K);
M.W    = zeros(cfg.mti.maxLength, K);
for i = 1:K
    B = cpf_band(cfg, prfHz(i), lambda);
    [w, M.caDb(i), edge] = cpf_design_filter(B, design.sigmaHz / prfHz(i), design.cnrDb, cfg);
    M.W(:, i) = w;
    M.G(i, :) = abs(B.S * w).^2;
    if strcmpi(cfg.prf.blindMode, 'notch')
        M.blind(i, :) = M.G(i, :) < 10^(cfg.prf.blindDb / 10);
    else
        M.blind(i, :) = B.dn.' <= edge + 0.5 / cfg.prf.fftLength;
    end
end
end
