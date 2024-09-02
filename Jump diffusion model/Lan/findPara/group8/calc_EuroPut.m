function y = calc_EuroPut(S0, T, r, sigma, K)
d1 = (log(S0 / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T));
d2 = d1 - sigma * sqrt(T);
y = K * exp(-r * T) .* (1 - normcdf(d2, 0, 1)) - S0 .* (1 - normcdf(d1, 0, 1));