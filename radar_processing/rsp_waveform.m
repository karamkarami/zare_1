function s = rsp_waveform(pulse, fs)
%RSP_WAVEFORM Complex baseband replica of one transmitted pulse.
%
%   s = rsp_waveform(pulse, fs)
%
%   pulse : one entry of P.pulse (see rsp_default_params)
%   fs    : sampling frequency [MHz]
%   s     : column vector, unit amplitude, length round(widthUs*fs)
%
%   Priority: pulse.samples (user replica) > pulse.code (phase code) > LFM.
%   An unmodulated pulse is an LFM with bwMHz = 0.

if isfield(pulse, 'samples') && ~isempty(pulse.samples)
    s = double(pulse.samples(:));
    return
end

L = round(pulse.widthUs * fs);
if L < 1
    error('rsp_waveform:length', 'Pulse of %.3g us is shorter than one sample at %.3g MHz.', ...
        pulse.widthUs, fs);
end
t = ((0:L-1)' - (L-1)/2) / fs;          % time centred on the pulse [us]

fc = getField(pulse, 'fcMHz', 0);
bw = getField(pulse, 'bwMHz', 0);
sl = getField(pulse, 'slope', 1);
k  = sl * bw / pulse.widthUs;            % chirp rate [MHz/us]
s  = exp(1j*2*pi*(fc*t + 0.5*k*t.^2));

code = getField(pulse, 'code', []);
if ~isempty(code)
    chip = floor((0:L-1)' * numel(code) / L) + 1;   % sample -> chip index
    c    = double(code(:));
    s    = s .* c(chip) ./ abs(c(chip));
end
end

function v = getField(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
