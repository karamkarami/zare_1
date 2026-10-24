function S = cpf_steering(fHz, prfHz, N)
%CPF_STEERING Slow-time steering vectors, one row per Doppler frequency.
%
%   S = cpf_steering(fHz, prfHz, N)    S(i, k+1) = exp(-j 2 pi f_i k / PRF)
%
%   With FIR taps w (y(m) = sum_k w(k+1) x(m-k)) the gain at f_i is |S(i,:) * w|^2.

k = 0:N-1;
S = exp(-2j * pi * fHz(:) * k / prfHz);
end
