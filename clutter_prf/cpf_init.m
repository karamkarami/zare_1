function st = cpf_init(cfg, prfHz)
%CPF_INIT State of the clutter-map MTI: empty map, first PRF set.
%
%   st = cpf_init(cfg)          the first PRF set is selected for cfg.clutter
%                               (no map yet), cfg.radar.rpm, stationary platform
%   st = cpf_init(cfg, prfHz)   start with this PRF set instead (no selection)
%   Then per block:  [y, st, info] = cpf_process(st, blk)
%   and per scan:    [prf, st, sel] = cpf_next_prf(st)

st.cfg    = cfg;
st.lambda = 299792458 / cfg.radar.fcHz;
nA = round(360 / cfg.map.azCellDeg);
nR = ceil(cfg.map.maxRangeM / cfg.map.rangeCellM);
st.map = cpf_map_empty(cfg, nA, nR);
st.acc = accEmpty(nA * nR, cfg.map.maxLag);
st.accPrf   = [];
st.noise    = cfg.noise.power;
st.omegaDps = [];                       % measured antenna speed [deg/s]
st.pos      = [0 0];                    % platform position at the next pulse
st.vel      = [0 0];
st.phi      = 0;                        % platform Doppler phase of the last pulse
st.tail     = [];                       % last pulses of the previous block
st.lastAz   = [];
st.scan     = 0;                        % revolutions seen in the data
st.banks    = {};
st.bankPrf  = [];

st.schedScan = 0;                       % scans scheduled by cpf_next_prf
st.selScan   = 0;
st.selSpeed  = 0;
if nargin < 2 || isempty(prfHz)
    st.sel = cpf_select_prfs(cfg, cpf_design_clutter(st, 0), st.lambda);
else
    st.sel = cpf_prf_eval(cfg, cpf_design_clutter(st, 0), sort(prfHz(:).'), st.lambda);
end
st.prfSet    = st.sel.prfHz;
st.prfNext   = 1;
end

function a = accEmpty(n, K)
a.s0 = zeros(n, 1);  a.n0 = zeros(n, 1);       % lag 0
a.s  = complex(zeros(n, K));  a.n = zeros(n, K);   % lags 1..K
a.m2 = zeros(n, 1);                            % platform spread^2, summed
end
