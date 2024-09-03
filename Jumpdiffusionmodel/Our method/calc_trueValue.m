function y = calc_trueValue(S0, D, n, tau, T, r, mu, sigma, K, alpha, J_lambda, J_m, J_theta)
CorMatrix = 0.3 * ones(D, D) + 0.7 * diag(ones(1, D));
temp = chol(CorMatrix);
CorL = temp';
V0 = calc_JumpEuroPut(S0, T, r, sigma, K, J_lambda, J_m, J_theta) * D;
J_r = mu - (sigma^2)/2 - (exp(J_m + (J_theta^2)/2) - 1) * J_lambda;
J_J = (rand([D, n]) >= exp(-J_lambda * tau));
sampleX = S0 * exp(J_r * tau * ones(D, n) ...
    + sigma * CorL * sqrt(tau) * normrnd(0, 1, D, n) ...
    + normrnd(0, 1, D, n) * J_theta .* J_J + J_m * J_J ...
    );
value = V0;
for j=1:D
    value = value - calc_JumpEuroPut(sampleX(j, :) , T - tau, r, sigma, K, J_lambda, J_m, J_theta);
end
value = sort(value, 'descend');
mm = ceil((1 - alpha) * n);
y = sum(value(1: mm), 'all') / mm;

