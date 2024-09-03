function y = calc_delta_basis(SS, TT, r, sigma, KK)
d1 = (log(SS / KK) + (r + sigma^2 / 2) * TT) / (sigma * sqrt(TT));
y = (normcdf(d1, 0, 1) - 1) .* SS;