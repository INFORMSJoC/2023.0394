replications = 50;

ave_gap = 0;
for repid = 1:replications
    filename = strcat("results\", num2str(repid), ".txt");
    fid = fopen(filename, 'r');
    res = fscanf(fid, "%f");
    ll = res(1);
    rr = res(2);
    ave_gap = ave_gap + (rr - ll);
    fclose(fid);
end
ave_gap = ave_gap / replications;
filename = strcat("ave_CI.txt");
fid = fopen(filename, 'w');
fprintf(fid, "%.6f\n", ave_gap);
fclose(fid);