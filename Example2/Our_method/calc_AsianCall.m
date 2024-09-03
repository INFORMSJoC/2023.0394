function y = calc_AsianCall(S0, prelogsum, delta_t, m, l1, r, sigma, K)
size_S0 = size(S0);
size_K = size(K);
new_S0 = repmat(S0, [1, 1, size_K(2)]);
new_prelogsum = repmat(prelogsum, [1, 1, size_K(2)]);
new_K = reshape(K, [1, size_K]);
new_K = repmat(new_K, [size_S0, 1]);

sigma2_ = sigma * sigma * delta_t;
tt = 1.0 * (m - l1) * (m - l1) / m / m;
tt = tt + (m - l1 - 1) * (m - l1) * (2 * m - 2 * l1 - 1) / 6.0 / (m^2);
sigma2_ = sigma2_ * tt;
mu2 = (r - sigma * sigma / 2) * delta_t;
mu2 = mu2 * ((m - l1) / m + (m - l1 - 1) * (m - l1) / 2.0 / m);
new_S0 = exp(new_prelogsum / m) .* power(new_S0, (m - l1) / m);
d1 = (log(new_S0 ./ new_K) + mu2 + sigma2_) / sqrt(sigma2_);
d2 = d1 - sqrt(sigma2_);
y =(new_S0 * exp(mu2 + sigma2_ / 2) .* normcdf(d1, 0, 1) - new_K .* normcdf(d2, 0, 1)) / exp(r * (m - l1) * delta_t);