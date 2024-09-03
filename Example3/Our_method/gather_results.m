sample_size_set = [1 5 10 20 40 50 60 80 100] * 10^4;
replications = 1000;
filename = strcat("results\true_value.txt");
fid = fopen(filename, 'r');
true_value = fscanf(fid, "%f");
fclose(fid);
res = zeros(replications, 7);
coverage_probs = zeros(length(sample_size_set), 1);
CI_ratios = zeros(length(sample_size_set), 1);

Our_point_MSEs = zeros(length(sample_size_set), 1);
Our_point_RRMSEs = zeros(length(sample_size_set), 1);
LS_point_MSEs = zeros(length(sample_size_set), 1);
LS_point_RRMSEs = zeros(length(sample_size_set), 1);
UB_outer = zeros(length(sample_size_set), 1);
LB_outer = zeros(length(sample_size_set), 1);
for ss = 1:length(sample_size_set)
    n = sample_size_set(ss);
    coverage = 0;
    ave_gap = 0;
    Our_point_MSE = 0;
    LS_point_MSE = 0;
    for repid = 1:replications
        filename = strcat("results\", num2str(n), "\", num2str(repid), ".txt");
        fid = fopen(filename, 'r');
        res(repid, :) = fscanf(fid, "%f");
        ll = res(repid, 6);
        rr = res(repid, 7);
        ave_gap = ave_gap + (rr - ll);
        coverage = coverage + (ll <= true_value) * (true_value <= rr);
		Our_point_MSE = Our_point_MSE + ((rr + ll) / 2 - true_value)^2;
		LS_point_MSE = LS_point_MSE + (res(repid, 5) - true_value)^2;
	if ll > true_value
            LB_outer(ss) = LB_outer(ss) + 1;
        end
        if rr < true_value
            UB_outer(ss) = UB_outer(ss) + 1;
        end
        fclose(fid);
    end
    ave_gap = ave_gap / replications;
	Our_point_MSE = Our_point_MSE / replications;
    Our_point_RRMSE = sqrt(Our_point_MSE) / true_value * 100;
	LS_point_MSE = LS_point_MSE / replications;
    LS_point_RRMSE = sqrt(LS_point_MSE) / true_value * 100;
    
    filename = strcat("results", num2str(n), ".csv");
    fid = fopen(filename, 'w');
    fprintf(fid, "%s,%.6f\n", "true_value", true_value);
    fprintf(fid, "%s, %.3f\n", "coverage_prob", coverage / replications);
    fprintf(fid, "%s, %.4f\n", "average LB", sum(res(:, 1), 'all') / replications);
    fprintf(fid, "%s, %.4f\n", "average UB", sum(res(:, 3), 'all') / replications);
    fprintf(fid, "%s, %.4f\n", "average LS", sum(res(:, 5), 'all') / replications);
    fprintf(fid, "%s, %.4f\n", "average CI_left", sum(res(:, 6), 'all') / replications);
    fprintf(fid, "%s, %.4f\n", "average CI_right", sum(res(:, 7), 'all') / replications);
    fprintf(fid, "%s, %.4f, %s, %.4f\n", "average width", ave_gap, "ratio to true value", ave_gap / true_value);
	fprintf(fid, "%s, %.6f, %s, %.4f\n", "our MSE of point estimation", Our_point_MSE, "our RRMSE of point estimation", Our_point_RRMSE);
	fprintf(fid, "%s, %.6f, %s, %.4f\n", "LS MSE of point estimation", LS_point_MSE, "LS RRMSE of point estimation", LS_point_RRMSE);
    fprintf(fid, "%s,%s,%s,%s,%s,%s,%s,%s\n", "LB estimate", "LB sigma2", "UB estimate", "UB sigma2", "LS_estimate_", "CI_left", "CI_right", "in CI");
    for repid = 1:replications
        for k=1:7
            fprintf(fid, "%.4f,", res(repid, k));
        end
        fprintf(fid, "%d\n", (res(repid, 6) <= true_value) * (true_value <= res(repid, 7)));
    end 
    fclose(fid);
    coverage_probs(ss) = coverage / replications;
    CI_ratios(ss) = ave_gap / true_value;
	
	Our_point_MSEs(ss) = Our_point_MSE;
    Our_point_RRMSEs(ss) = Our_point_RRMSE;
	LS_point_MSEs(ss) = LS_point_MSE;
    LS_point_RRMSEs(ss) = LS_point_RRMSE;
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
    fprintf(fid, "%.8f ", Our_point_MSEs(ss));
end
fprintf(fid, '\n');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%.4f ", Our_point_RRMSEs(ss));
end
fprintf(fid, '\n');

for ss = 1:length(sample_size_set)
    fprintf(fid, "%.8f ", LS_point_MSEs(ss));
end
fprintf(fid, '\n');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%.4f ", LS_point_RRMSEs(ss));
end
fprintf(fid, '\n');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%d ", LB_outer(ss));
end
fprintf(fid, '\n');
for ss = 1:length(sample_size_set)
    fprintf(fid, "%d ", UB_outer(ss));
end
fprintf(fid, '\n');
fclose(fid);