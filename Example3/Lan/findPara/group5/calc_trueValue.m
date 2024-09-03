function y = calc_trueValue(S0, D, n, m, l1, T, r, mu, sigma, K, H, alpha)
CorMatrix = 0.3 * ones(D, D) + 0.7 * diag(ones(1, D));
temp = chol(CorMatrix);
CorL = temp';
V0 = calc_V0(S0, T, r, sigma, K, H) * D;
delta_t = T / m;
tau = l1 * delta_t;
samplePath = zeros(D, n, l1);
samplePath(:, :, 1) = exp((mu - sigma^2 / 2) * delta_t ...
    + sigma * sqrt(delta_t) * CorL * normrnd(0, 1, D, n)) * S0;
for l = 2:l1
    samplePath(:, :, l) = exp((mu - sigma^2 / 2) * delta_t ...
    + sigma * sqrt(delta_t) * CorL * normrnd(0, 1, D, n)) .* samplePath(:, :, l - 1);
end
sampleMin = zeros(D, n, l1);
UnifRand = rand(D, n, l1);
sampleMin(:, :, 1) = exp((log(S0 * samplePath(:, :, 1)) - sqrt(power(log(samplePath(:, :, 1) / S0), 2) - 2 * sigma^2 * delta_t * log(UnifRand(:, :, 1))))/2);
for l = 2:l1
    sampleMin(:, :, l) = exp((log(samplePath(:, :, l-1) .* samplePath(:, :, l)) - sqrt(power(log(samplePath(:, :, l) ./ samplePath(:, :, l-1)), 2) - 2 * sigma^2 * delta_t * log(UnifRand(:, :, l))))/2);
end
sampleMin = min(sampleMin(:, :, :), [], 3);
value = calc_BarrierCall(samplePath(:, :, l1), T - tau, r, sigma, K(1), H) .* (sampleMin(:, :) >= H);
value = calc_BarrierCall(samplePath(:, :, l1), T - tau, r, sigma, K(2), H) .* (sampleMin(:, :) >= H) + value;
value = calc_BarrierCall(samplePath(:, :, l1), T - tau, r, sigma, K(3), H) .* (sampleMin(:, :) >= H) + value;
clear samplePath sampleMin;
value = V0 - sum(value, 1);
value = sort(value, 'descend');
sz = ceil((1 - alpha) * n);
y = sum(value(1: sz)) / sz;
% 28.5867
