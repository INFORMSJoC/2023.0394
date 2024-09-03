function y = calc_V0(S0, T, r, sigma, K, H)
y = sum(calc_BarrierCall(S0, T, r, sigma, K, H), 'all');