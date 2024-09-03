function y = calc_JumpEuroPut(S0, T, r, sigma, K, J_lambda, J_m, J_theta)
ret_y = 0;
facs = 1;
for k=0:20
    if (k >= 2)
        facs = facs * k;
    end
    J_sigma = sqrt((sigma^2) + k * (J_theta^2) / T);
    J_S0 = S0 * exp(k * (J_m + (J_theta^2)/2) - J_lambda * (exp(J_m + (J_theta^2)/2) - 1) * T);
    ret_y = ret_y + ((J_lambda * T)^k) * (calc_EuroPut(J_S0, T, r, J_sigma, K)) / facs;
end
y = exp(-J_lambda * T) * ret_y;