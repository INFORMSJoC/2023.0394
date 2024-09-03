function y = calc_trueValue(S0, D, n, tau, T, r, mu, sigma, K, alpha)
CorMatrix = 0.3 * ones(D, D) + 0.7 * diag(ones(1, D));
temp = chol(CorMatrix);
CorL = temp';
V0 = calc_V0(S0, T, r, sigma, K) * D;
sampleX = exp((mu - sigma^2 / 2) * tau * ones(D, n) ...
    + sigma * sqrt(tau) * CorL * normrnd(0, 1, D, n)) * S0;
sampleX = sampleX';
value = calc_EuroCall(sampleX, T - tau, r, sigma, K(1));
value = value + calc_EuroCall(sampleX, T - tau, r, sigma, K(2));
value = value + calc_EuroCall(sampleX, T - tau, r, sigma, K(3));
clear sampleX;
value = V0 - sum(value, 2);
value = sort(value, 'descend');
m = ceil((1 - alpha) * n);
y = sum(value(1: m)) / m;

% 1.498144795195251e+02
