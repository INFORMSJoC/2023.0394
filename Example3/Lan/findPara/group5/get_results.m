replications = 20;
p = 0.05;
alpha = 0.05;
D = 3;
T = 1;
m = 50;
H = 80;
l1 = 3;
S0 = 100;
K = [90 100 110];
mu = 8/100;
r = 5/100;
sigma = 0.15;

alphao = alpha / 2; % alpha outer
alphai = alpha - alphao;
alphas = alphai * 0.4;
alphalo = alphai * 0.3;
alphahi = alphai * 0.3;
c = exp(-(1/2)*chi2inv(1-alpha,1));


Budget = 10^6;
n0 = fix(Budget ^ (3/5))
m0 = fix(Budget * 0.4 / n0)

L = zeros(1,n0);
for l = 1:n0
    if (n0*log(n0)+l*(log(p)-log(l)) +(n0-l)*(log(1-p)-log(n0-l)) >=log(c) )
        L(l) = l;
    end
end
lmin = min(L(L~=0));
lmax = max(L);

clear L
delta = zeros(1,n0);
for l = lmin:lmax
    delta(l) = sqrt(MAX(l,p,c,n0));
end
results = zeros(replications, 2);
for repid = 1:replications
    filename = strcat("results\", num2str(repid), ".txt");
    fid = fopen(filename, 'w');
    results(repid,:) = CIofCVaR(D, S0, K, H, mu, r, T, sigma, l1, m, alphas,alphahi,alphalo,p,Budget,c,n0,m0,delta,lmin,lmax);
    fprintf(fid,"%.4f %.4f\n", results(repid,1), results(repid,2));
    fclose(fid);
end
