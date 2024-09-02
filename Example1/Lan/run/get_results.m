sample_size_set = [5 10 20 40 50 60 80 100] * 10^4;
replications = 1000;
p = 0.05;
alpha = 0.05;
D = 20;
T = 1;
m = 50;
h = T/m;
l1 = 3;
tau = l1*h;
S0 = 100;
K = [90 100 110];
mu = 8/100;
r = 5/100;
sigma = 0.15;

alphao = alpha / 2;
alphai = alpha - alphao;
alphas = alphai * 0.4;
alphalo = alphai * 0.3;
alphahi = alphai * 0.3;
c = exp(-(1/2)*chi2inv(1-alpha,1));

for ss = 1:length(sample_size_set)
    Budget = sample_size_set(ss)

    n0 = fix(Budget^(1/2));
    m0 = fix(Budget * 0.4 / n0);
	if n0 < 800
		n0 = 800;
	end
	if m0 < 30
		m0 = 30;
    end
    n0
    m0
	if n0 * m0 >= Budget
		continue;
    end

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
        filename = strcat("results\", num2str(Budget), "\", num2str(repid), ".txt");
        fid = fopen(filename, 'w');
        results(repid,:) = CIofCVaR(D, S0, K, mu, r, T, sigma, tau, alpha, alphao, alphai,alphas,alphahi,alphalo,p,Budget,c,n0,m0,delta,lmin,lmax);
        fprintf(fid,"%.4f %.4f\n", results(repid,1), results(repid,2));
        fclose(fid);
    end
end

