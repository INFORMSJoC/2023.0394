%sample_size_set = [1 2 4 6 8 10] * 10^5;
%sample_size_set = [10000 50000 sample_size_set];
sample_size_set = [400000]
S0 = 100;
D = 20;
T = 1;
m = 50;
l1 = 3; %tau = l1/m*T
r = 0.05;
mu = 0.08;
sigma = 0.15;
K = [90 100 110];
alpha = 0.95;
z_beta = 0.05;
z_coef = norminv(1 - z_beta / 2);
replications = 1000;
need_true_value = 0;
results = zeros(replications, 7);
% save true value
if need_true_value > 0
    filename = strcat("results\true_value.txt");
    fid = fopen(filename, 'w');
    true_value = calc_trueValue(S0, D, 10^7, l1 / m * T, T, r, mu, sigma, K, alpha);
    fprintf(fid, "%.6f\n", true_value);
    fclose(fid);
end

for ss = 1:length(sample_size_set)
    n = sample_size_set(ss);
    for repid = 882:replications
        filename = strcat("results\", num2str(n), "\", num2str(repid), ".txt");
        fid = fopen(filename, 'w');
        results(repid, 1:5) = construct_CI(S0, D, n, T, m, l1, r, mu, sigma, K, alpha);
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
