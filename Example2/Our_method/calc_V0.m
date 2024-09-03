function y = calc_V0(S0, delta_t, m, r, sigma, K)
y = sum(calc_AsianCall(S0, zeros(size(S0)), delta_t, m, 0, r, sigma, K), 'all');