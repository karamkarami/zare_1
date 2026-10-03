function [s, spc] = rsp_waveform(pulse, fs)
%RSP_WAVEFORM Complex baseband samples of one transmitted pulse.
%
%   [s, spc] = rsp_waveform(pulse, fs)
%
%   pulse : one entry of P.pulse (see rsp_default_params)
%   fs    : sampling frequency [MHz]
%   s     : column vector, unit amplitude
%   spc   : samples per chip (code) or samples per 1/bandwidth (LFM)
%
%   pulse.type
%     'lfm'  : widthUs long, bwMHz swept, slope +1/-1 (bwMHz = 0: plain pulse)
%     'code' : pulse.code holds the chip values (any length, binary or
%              polyphase). With codeRate = 'chip' every chip lasts 1/bwMHz us
%              (bwMHz = 1 -> 1 us chips); with codeRate = 'sample' the vector
%              is already sampled at fs.
%   pulse.fcMHz shifts the pulse inside the baseband.
%   pulse.samples (not empty) overrides everything: user samples at fs.

if isfield(pulse, 'samples') && ~isempty(pulse.samples)
    s   = double(pulse.samples(:));
    spc = 1;
    return
end

bw = getField(pulse, 'bwMHz', 0);
switch lower(getField(pulse, 'type', 'lfm'))
    case 'lfm'
        L = round(pulse.widthUs * fs);
        if L < 1
            error('rsp_waveform:length', 'Pulse shorter than one sample.');
        end
        t = ((0:L-1)' - (L-1)/2) / fs;                % centred time [us]
        k = getField(pulse, 'slope', 1) * bw / pulse.widthUs;
        s = exp(1j*pi*k*t.^2);
        spc = fs / max(bw, 1/pulse.widthUs);
    case 'code'
        c = double(pulse.code(:));
        if isempty(c)
            error('rsp_waveform:code', 'pulse.code is empty.');
        end
        c = c ./ abs(c);
        if strcmpi(getField(pulse, 'codeRate', 'chip'), 'sample')
            s   = c;
            spc = 1;
        else
            if bw <= 0
                error('rsp_waveform:bw', 'bwMHz (chip rate) must be > 0 for a code.');
            end
            spc  = fs / bw;                            % samples per chip
            L    = round(numel(c) * spc);
            chip = floor((0:L-1)' / spc) + 1;          % sample -> chip
            s    = c(min(chip, numel(c)));
        end
    otherwise
        error('rsp_waveform:type', 'Unknown pulse type "%s".', pulse.type);
end

fc = getField(pulse, 'fcMHz', 0);
if fc ~= 0
    s = s .* exp(1j*2*pi*fc*(0:numel(s)-1)'/fs);
end
end

function v = getField(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
