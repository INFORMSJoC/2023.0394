sample_size_set = [1 2 4 8] * 10^5;
replications = 100;
filename = strcat("results\true_value.txt");
fid = fopen(filename, 'r');
true_value = fscanf(fid, "%f");
fclose(fid);

coverage_probs = zeros(length(sample_size_set), 1);
CI_ratios = zeros(length(sample_size_set), 1);
point_MSEs = zeros(length(sample_size_set), 1);
point_RRMSEs = zeros(length(sample_size_set), 1);

for ss = 1:length(sample_size_set)
    budget = sample_size_set(ss);
    res = zeros(replications, 2);
    coverage = 0;
    ave_gap = 0;
    point_MSE = 0;
    for repid = 1:replications
        filename = strcat("results\", num2str(budget), "\", num2str(repid), ".txt");
        fid = fopen(filename, 'r');
        res(repid, :) = fscanf(fid, "%f");
        ll = res(repid, 1);
        rr = res(repid, 2);
        ave_gap = ave_gap + (rr - ll);
        coverage = coverage + (ll <= true_value) * (true_value <= rr);
        point_MSE = point_MSE + ((rr + ll) / 2 - true_value)^2;
        fclose(fid);
    end
    ave_gap = ave_gap / replications;
    point_MSE = point_MSE / replications;
    point_RRMSE = sqrt(point_MSE) / true_value * 100;

    filename = strcat("results", num2str(budget), ".csv");
    fid = fopen(filename, 'w');
    fprintf(fid, "%s,%.6f\n", "true_value", true_value);
    fprintf(fid, "%s, %.3f\n", "coverage_prob", coverage / replications);
    fprintf(fid, "%s, %.4f, %s, %.4f\n", "average width", ave_gap, "ratio to true value", ave_gap / true_value);
    fprintf(fid, "%s, %.6f, %s, %.4f\n", "MSE of point estimation", point_MSE, "RRMSE of point estimation", point_RRMSE);
    fprintf(fid, "%s, %.4f\n", "average CI_left", sum(res(:, 1)) / replications);
    fprintf(fid, "%s, %.4f\n", "average CI_right", sum(res(:, 2)) / replications);
    fprintf(fid, "%s,%s,%s\n", "CI_left", "CI_right", "in CI");
    for repid = 1:replications
        for k=1:2
            fprintf(fid, "%.4f,", res(repid, k));
        end
        fprintf(fid, "%d\n", (res(repid, 1) <= true_value) * (true_value <= res(repid, 2)));
    end 
    fclose(fid);
    
    coverage_probs(ss) = coverage / replications;
    CI_ratios(ss) = ave_gap / true_value;
    point_MSEs(ss) = point_MSE;
    point_RRMSEs(ss) = point_RRMSE;
end
fid = fopen("summary.txt", 'w');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%d ", sample_size_set(ss));
end
fprintf(fid, '\n');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%.3f ", coverage_probs(ss));
end
fprintf(fid, '\n');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%.4f ", CI_ratios(ss));
end
fprintf(fid, '\n');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%.8f ", point_MSEs(ss));
end
fprintf(fid, '\n');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%.4f ", point_RRMSEs(ss));
end
fprintf(fid, '\n');

fclose(fid);