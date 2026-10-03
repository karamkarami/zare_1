function d = rsp_code_mmf(code, len)
%RSP_CODE_MMF Least-squares mismatched filter (decoder) of a phase code.
%
%   d = rsp_code_mmf(code, len)
%
%   code : chip values of the code
%   len  : decoder length in chips (>= numel(code)); longer = lower sidelobes
%   d    : FIR taps at chip rate (convolution form, P.pulse(k).decoderForm = 'fir'),
%          scaled so the peak response is numel(code)
%
%   The taps minimise the energy of conv(code, d) away from the main peak
%   (integrated sidelobe level). Example for Barker 13:
%       P.pulse(1).decoder = rsp_code_mmf(rsp_code('barker', 13), 39);

code = code(:);
N    = numel(code);
if len < N
    error('rsp_code_mmf:len', 'The decoder must be at least as long as the code.');
end
L = N + len - 1;                                    % output length
C = zeros(L, len);                                  % convolution matrix
for k = 1:len
    C(k:k+N-1, k) = code;
end
e = zeros(L, 1);
e((L + 1)/2 + mod(L + 1, 2)/2) = N;                 % desired: peak in the middle
d = (C \ e).';
end
