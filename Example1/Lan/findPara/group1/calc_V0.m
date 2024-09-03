function y = calc_V0(S0, T, r, sigma, K)
y = sum(calc_EuroCall(S0, T, r, sigma, K), 'all');