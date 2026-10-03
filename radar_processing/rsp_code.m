function c = rsp_code(name, N)
%RSP_CODE Common phase codes (row vector, unit magnitude).
%
%   c = rsp_code(name, N)
%
%   'barker' : N = 2, 3, 4, 5, 7, 11, 13
%   'mls'    : maximal-length binary sequence, N = 2^m - 1 (m = 3..12)
%   'frank'  : N = M^2
%   'p3'     : any N
%   'p4'     : any N
%
%   Use your own codes directly in P.pulse(k).code; this is only a helper.

switch lower(name)
    case 'barker'
        switch N
            case 2,  c = [1 -1];
            case 3,  c = [1 1 -1];
            case 4,  c = [1 1 -1 1];
            case 5,  c = [1 1 1 -1 1];
            case 7,  c = [1 1 1 -1 -1 1 -1];
            case 11, c = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
            case 13, c = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
            otherwise
                error('rsp_code:barker', 'No Barker code of length %d.', N);
        end
    case 'mls'
        m = round(log2(N + 1));
        if 2^m - 1 ~= N || m < 3 || m > 12
            error('rsp_code:mls', 'MLS length must be 2^m - 1 with m = 3..12.');
        end
        taps = {[3 2], [4 3], [5 3], [6 5], [7 6], [8 6 5 4], [9 5], [10 7], [11 9], ...
                [12 11 10 4]};
        t   = taps{m - 2};
        reg = ones(1, m);
        c   = zeros(1, N);
        for i = 1:N
            c(i) = reg(m);
            fb   = mod(sum(reg(t)), 2);
            reg  = [fb reg(1:m-1)];
        end
        c = 1 - 2*c;                                     % {0,1} -> {+1,-1}
    case 'frank'
        M = round(sqrt(N));
        if M^2 ~= N
            error('rsp_code:frank', 'Frank code length must be a square.');
        end
        [i, k] = ndgrid(0:M-1, 0:M-1);
        ph = 2*pi/M * (i .* k);
        c  = exp(1j * ph(:)).';
    case 'p3'
        n = 0:N-1;
        c = exp(1j*pi*n.^2/N);
    case 'p4'
        n = 0:N-1;
        c = exp(1j*(pi*n.^2/N - pi*n));
    otherwise
        error('rsp_code:name', 'Unknown code "%s".', name);
end
end
