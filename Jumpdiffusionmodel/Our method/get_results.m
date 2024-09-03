sample_size_set = [1 2 4 8] * 10^5;
S0 = 100;
D = 20;
T = 1/4;
m = 50;
l1 = 3;
r = 0.05;
mu = 0.08;
sigma = 0.20;
J_lambda = 1;
J_m = 0.03;
J_theta = 0.15;
K = 120;
alpha = 0.90;
z_beta = 0.05;
z_coef = norminv(1 - z_beta / 2);
replications = 100;
need_true_value = 0;
results = zeros(replications, 7);
% save true value
if need_true_value > 0
    filename = strcat("results\true_value.txt");
    fid = fopen(filename, 'w');
    true_value = calc_trueValue(S0, D, 10^6, T/m*l1, T, r, mu, sigma, K, alpha, J_lambda, J_m, J_theta);
    fprintf(fid, "%.6f\n", true_value);
    fclose(fid);
end

for ss = 1:length(sample_size_set)
    n = sample_size_set(ss);
    for repid = 1:replication
        filename = strcat("results\", num2str(n), "\", num2str(repid), ".txt");
        fid = fopen(filename, 'w');
        results(repid, 1:5) = construct_CI(S0, D, n, T, m, l1, mu, r, sigma, K, J_lambda, J_m, J_theta, alpha);
        ll = results(repid, 1) - sqrt(results(repid, 2) / n) * z_coef;
        rr = results(repid, 3) + sqrt(results(repid, 4) / n) * z_coef;
        results(repid, 6) = ll;
        results(repid, 7) = rr;
        for k=1:7
            fprintf(fid, "%.4f ", results(repid, k));
        end
        fprintf(fid, "\n");
        fclose(fid);
    end
end
