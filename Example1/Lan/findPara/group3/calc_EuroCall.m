function y = calc_EuroCall(S0, T, r, sigma, K)
size_S0 = size(S0);
size_K = size(K);
new_S0 = repmat(S0, [1, 1, size_K(2)]);
new_K = reshape(K, [1, size_K]);
new_K = repmat(new_K, [size_S0, 1]);
d1 = (log(new_S0 ./ new_K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T));
d2 = d1 - sigma * sqrt(T);
y = new_S0 .* normcdf(d1, 0, 1) - new_K * exp(-r * T) .* normcdf(d2, 0, 1);
clear d1 d2 size_S0 size_K new_S0 new_K;