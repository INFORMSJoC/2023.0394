sample_size_set = [1 5 50] * 10^4;
S0 = 100;
D = 3;
T = 1;
m = 50;
l1 = 3; %tau = l1/m*T
r = 0.05;
mu = 0.08;
sigma = 0.15;
K = [90 100 110];
H = 80;
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
    true_value = calc_trueValue(S0, D, 10^7, m, l1, T, r, mu, sigma, K, H, alpha);
    fprintf(fid, "%.6f\n", true_value);
    fclose(fid);
end
for ss = 1:length(sample_size_set)
    n = sample_size_set(ss);
    for repid = 1:replications
        filename = strcat("results\", num2str(n), "\", num2str(repid), ".txt");
        fid = fopen(filename, 'w');
        results(repid, 1:5) = construct_CI(S0, D, n, m, l1, T, r, mu, sigma, K, H, alpha);
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
