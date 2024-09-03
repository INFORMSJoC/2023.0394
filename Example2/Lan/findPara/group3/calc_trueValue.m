function y = calc_trueValue(S0, D, n, l1, m, T, r, mu, sigma, K, alpha)
CorMatrix = 0.3 * ones(D, D) + 0.7 * diag(ones(1, D));
temp = chol(CorMatrix);
CorL = temp';
delta_t = T / m;
V0 = calc_V0(S0, delta_t, m, r, sigma, K) * D;

samplePath = zeros(D, n, l1);
samplePath(:, :, 1) = (mu - sigma^2 / 2) * delta_t ...
    + sigma * sqrt(delta_t) * CorL * normrnd(0, 1, D, n) + log(S0);
for l = 2:l1
    samplePath(:, :, l) = (mu - sigma^2 / 2) * delta_t ...
    + sigma * sqrt(delta_t) * CorL * normrnd(0, 1, D, n) + samplePath(:, :, l - 1);
end
prelogsum = zeros(D, n);
for l = 1:l1
    prelogsum = prelogsum + samplePath(:, :, l);
end
X = exp(samplePath(:, :, l1));
value = calc_AsianCall(X, prelogsum, delta_t, m, l1, r, sigma, K(1));
value = calc_AsianCall(X, prelogsum, delta_t, m, l1, r, sigma, K(2)) + value;
value = calc_AsianCall(X, prelogsum, delta_t, m, l1, r, sigma, K(3)) + value;
value = V0 - sum(value, 1);
value = sort(value, 'descend');
m = ceil((1 - alpha) * n);
y = sum(value(1: m)) / m;

% 23.8082
