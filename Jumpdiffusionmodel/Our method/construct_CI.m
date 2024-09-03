function y = construct_CI(S0, D, n, T, m, l1, mu, r, sigma, K, J_lambda, J_m, J_theta, alpha)
%lambda=1
il = [-1e10 -0.16223273 -0.09624318 -0.04866007 -0.00800206 0.03000000 0.06800207 0.10866008 0.15624319 0.22223274];
ir = [-0.16223273 -0.09624318 -0.04866007 -0.00800206 0.03000000 0.06800207 0.10866008 0.15624319 0.22223274 1e10];
yk = [-0.16223273 -0.09624318 -0.04866007 -0.00800206 0.03000000 0.03000000 0.06800207 0.10866008 0.15624319 0.22223274];

CorMatrix = 0.3 * ones(D, D) + 0.7 * diag(ones(1, D));
temp = chol(CorMatrix);
CorL = temp';
V0 = calc_JumpEuroPut(S0, T, r, sigma, K, J_lambda, J_m, J_theta);
delta_t = T / m;
N01 = normrnd(0, 1, D, n, m);
J_N01 = normrnd(0, 1, D, n, m);
J_J = (rand([D, n, m]) >= exp(-J_lambda * delta_t));
J_r1 = mu - (sigma^2)/2 - (exp(J_m + (J_theta^2)/2) - 1) * J_lambda;
J_r2 = r - (sigma^2)/2 - (exp(J_m + (J_theta^2)/2) - 1) * J_lambda;
samplePath = zeros(D, n, m);
samplePath0 = zeros(D, n, m);

samplePath0(:, :, 1) = S0 * exp(J_r1 * delta_t * ones(D, n)+ sigma * CorL * sqrt(delta_t) * N01(:, :, 1));
samplePath(:, :, 1) = samplePath0(:, :, 1) .* exp(J_N01(:, :, 1) * J_theta .* J_J(:, :, 1) + J_m * J_J(:, :, 1));
for l = 2:l1
    samplePath0(:, :, l) = samplePath(:, :, l - 1) .* exp(J_r1 * delta_t * ones(D, n)+ sigma * CorL * sqrt(delta_t) * N01(:, :, l));
    samplePath(:, :, l) = samplePath0(:, :, l) .* exp(J_N01(:, :, l) * J_theta .* J_J(:, :, l) + J_m * J_J(:, :, l));
end
for l = (l1 + 1):m
    samplePath0(:, :, l) = samplePath(:, :, l - 1) .* exp(J_r2 * delta_t * ones(D, n)+ sigma * CorL * sqrt(delta_t) * N01(:, :, l));
    samplePath(:, :, l) = samplePath0(:, :, l) .* exp(J_N01(:, :, l) * J_theta .* J_J(:, :, l) + J_m * J_J(:, :, l));
end
payoff = max(K - samplePath(:, :, m), 0);

% construct LB
sampleY = V0 * ones(D, n) - payoff * exp(-r * delta_t * (m - l1));
sampleY = sum(sampleY', 2);
varphi = zeros(n, D * 9 + 1);
varphi(:, 1) = 1;
idx = 1;
for j=1:D
    mD = ones(n, 1);
    mJ = ones(n, 1);
    for l = 1:l1
		tmD = exp(J_r1 * delta_t * ones(D, n) + sigma * CorL * sqrt(delta_t) * N01(:, :, l));
        mD = mD .* (tmD(j, :)');
		tmJ = exp(J_N01(:, :, l) * J_theta .* J_J(:, :, l) + J_m * J_J(:, :, l));
        mJ = mJ .* (tmJ(j, :)');
    end
    idx = idx + 1;
    varphi(:, idx) = mD;
    idx = idx + 1;
    varphi(:, idx) = mD .* mD;
	idx = idx + 1;
    varphi(:, idx) = mD .* mD .* mD;
    idx = idx + 1;
    varphi(:, idx) = mJ;
    idx = idx + 1;
    varphi(:, idx) = mJ .* mJ;
	idx = idx + 1;
    varphi(:, idx) = mJ .* mJ .* mJ;
    idx = idx + 1;
    varphi(:, idx) = mD .* mJ;
    idx = idx + 1;
    varphi(:, idx) = mD .* mD .* mJ .* mJ;
	idx = idx + 1;
    varphi(:, idx) = mD .* mD .* mD .* mJ .* mJ .* mJ;
end
beta = varphi \ sampleY; %(4D+1)*1
ybar = varphi * (beta);
aveY = sum(sampleY, 'all') / n;
sortedY = [ybar, sampleY]; % (n)*(2)
sortedY = sortrows(sortedY, 'ascend');

k_star = floor(n * alpha);
alpha_rev = 1.0 / (1 - alpha);
LS_estimate = sum(sortedY((k_star + 1): n, 1), 'all') / (n - k_star);

G_bar = (aveY - alpha * sortedY(1, 1)) * alpha_rev;
LB_u = sampleY(1, 1);
LB_estimate = G_bar;

for k=2:k_star
    G_bar = (((k - 1) / n - alpha) * sortedY(k, 1) - ((k - 2) / n - alpha) * sortedY(k - 1, 1)) ...
        * alpha_rev - alpha_rev * sortedY(k - 1, 2) / n + G_bar;
    if G_bar < LB_estimate
        LB_estimate = G_bar;
        LB_u = sortedY(k, 1);
    end
end
G_bar = G_bar + alpha_rev * (sortedY(k_star, 1) - sortedY(k_star, 2)) / n;
if G_bar < LB_estimate
    LB_estimate = G_bar;
    LB_u = sortedY(k_star, 1);
end
for k=(k_star+1):n
    G_bar = ((k / n - alpha) * sortedY(k, 1) - ((k - 1) / n - alpha) * sortedY(k - 1, 1)) ...
        * alpha_rev - alpha_rev * sortedY(k, 2) / n + G_bar;
    if G_bar < LB_estimate
        LB_estimate = G_bar;
        LB_u = sortedY(k, 1);
    end
end
clear sortedY;

% construct CLT for LB
Psi_n = (varphi') * varphi / n;
Phi_n = (sampleY') * varphi / n;
epsilon = std(ybar) * 1.06 * n^(-1/5);
Upsilon_n = (sampleY - LB_u) .* (1 - cos(ybar - LB_u)) .* (abs(ybar - LB_u) <= 2 * pi * epsilon) / (4 * pi * epsilon);
Upsilon_n = ((Upsilon_n') * varphi)' / n;
v2 = 0;
m1 = Phi_n / Psi_n;
m2 = Psi_n \ Upsilon_n;

for k = 1:n
    Psi_k = varphi(k, :);
    Psi_k = (Psi_k') * (Psi_k);
    Phi_k = sampleY(k, 1) * varphi(k, :);
    v = -m1 * Psi_k * m2;
    v = v + Phi_k * m2;
    v = v + (sampleY(k, 1) - LB_u) * (ybar(k, 1) >= LB_u);
    v2 = v2 + v * v;
end
LB_sigma2 = v2 / n - (sum((sampleY - LB_u) .* (ybar >= LB_u), 'all') / n) ^ 2;
LB_sigma2 = LB_sigma2 / ((1 - alpha)^2);


% contruct UB
for l = l1:(m - 1)
    coBM = sqrt(delta_t) * CorL * N01(:, :, l + 1);
    N01(:, :, l + 1) = coBM;
end

rr = 0.1 * (1 - exp(-J_lambda * delta_t)) - 0.01 * (1 - exp(-J_lambda * delta_t))^2;
delta_martingale = zeros(n, 1);
for l = l1:(m - 1)
    sampleY = V0 * ones(D, n) - payoff * exp(-r * delta_t * (m - l - 1));
    sampleY = sampleY';
    for j = 1:D
        if (l < m - 1)
            varphi = ones(n, 3);
            varphi(:, 1) = calc_delta_basis(samplePath0(j, :, l)', delta_t, r, sigma, K);
            varphi(:, 2) = calc_delta_basis(samplePath0(j, :, l)', delta_t * (m - l), r, sigma, K);
        else
            varphi = ones(n, 2);
            varphi(:, 1) = calc_delta_basis(samplePath0(j, :, l)', delta_t, r, sigma, K);
        end
        beta = varphi \ sampleY(:, j);
        residual_y = sampleY(:, j) - varphi * beta;
        coBM = N01(j, :, l + 1)';
        beta = varphi \ (residual_y .* coBM / delta_t);
        H = varphi * beta;
        delta_martingale = delta_martingale + H .* coBM;

        J_N = J_N01(j, :, l + 1) * J_theta + J_m;
        c1 = calc_EuroPut(samplePath0(j, :, l)', delta_t, r, sigma, K);
        if (l < m - 1)
            c2 = calc_EuroPut(samplePath0(j, :, l)', delta_t * (m - l), r, sigma, K);
        end 
        for ii = 1:10
            in_interval = (il(ii) <= J_N) .* (J_N <= ir(ii));
            coeff = (J_J(j, :, l + 1) .* in_interval) - (1 - exp(-J_lambda * delta_t)) * 0.1;
            coeff = coeff';
            if (l < m - 1)
                varphi = ones(n, 3);
                varphi(:, 1) = calc_EuroPut(exp(yk(ii)) * samplePath0(j, :, l)', delta_t, r, sigma, K) - c1;
                varphi(:, 2) = calc_EuroPut(exp(yk(ii)) * samplePath0(j, :, l)', delta_t * (m - l), r, sigma, K) - c2;
            else
                varphi = ones(n, 2);
                varphi(:, 1) = calc_EuroPut(exp(yk(ii)) * samplePath0(j, :, l)', delta_t, r, sigma, K) - c1;
            end
            beta = varphi \ sampleY(:, j);
            residual_y = sampleY(:, j) - varphi * beta;
            beta = varphi \ (residual_y .* coeff / rr);
            HH = varphi * beta;
            delta_martingale = delta_martingale + HH .* coeff;
        end
    end
end
sampleY = V0 * ones(D, n) - payoff * exp(-r * delta_t * (m - l1));
sampleY = (sum(sampleY, 1))';
z = sampleY - delta_martingale;
z = sort(z, 'descend');
sz = ceil((1 - alpha) * n);
UB_estimate = sum(z(1 : sz), 'all') / sz;

% construct CLT for UB
UB_u = z(sz);
z = sampleY - delta_martingale - UB_u;
signs = (z >= 0);
v2 = z .* signs;

for l = l1:(m - 1)
    sampleY = V0 * ones(D, n) - payoff * exp(-r * delta_t * (m - l - 1));
    sampleY = sampleY';
    for j = 1:D
        if (l < m - 1)
            varphi = ones(n, 3);
            varphi(:, 1) = calc_delta_basis(samplePath0(j, :, l)', delta_t, r, sigma, K);
            varphi(:, 2) = calc_delta_basis(samplePath0(j, :, l)', delta_t * (m - l), r, sigma, K);
            dd = 3;
        else
            varphi = ones(n, 2);
            varphi(:, 1) = calc_delta_basis(samplePath0(j, :, l)', delta_t, r, sigma, K);
            dd = 2;
        end
        coBM = N01(j, :, l + 1)';
        Psi_n = (varphi') * varphi / n;
        %yy = sampleY(:, j) .* coBM / delta_t;
        yy = sampleY(:, j);
        Phi_n = (yy') * varphi / n;
        Upsilon_n = (varphi') * (coBM .* signs) / n;
        ml = Phi_n / Psi_n;
        mr = Psi_n \ Upsilon_n;
        
        mat1 = repmat(reshape(varphi, n, dd, 1), [1, 1, dd]);
        mat2 = repmat(reshape(varphi, n, 1, dd), [1, dd, 1]);
        mat1 = mat1 .* mat2;
        mat1 = permute(mat1, [2 3 1]);
        mat1 = reshape(mat1, dd, dd * n);
        mat1 = ml * mat1;
        mat1 = reshape(mat1, dd, n);
        mat1 = (mat1') * mr;
        v2 = v2 + mat1;
        mat1 = repmat(yy, [1, dd]);
        mat1 = mat1 .* varphi;
        mat1 = mat1 * mr;
        v2 = v2 - mat1;        
    end
end


for l = l1:(m - 1)
    sampleY = V0 * ones(D, n) - payoff * exp(-r * delta_t * (m - l - 1));
    sampleY = sampleY';
    for j = 1:D
        J_N = J_N01(j, :, l + 1) * J_theta + J_m;
        c1 = calc_EuroPut(samplePath0(j, :, l)', delta_t, r, sigma, K);
        if (l < m - 1)
            c2 = calc_EuroPut(samplePath0(j, :, l)', delta_t * (m - l), r, sigma, K);
        end 
        for ii = 1:10
            in_interval = (il(ii) <= J_N) .* (J_N <= ir(ii));
            coeff = (J_J(j, :, l + 1) .* in_interval) - (1 - exp(-J_lambda * delta_t)) * 0.1;
            coeff = coeff';
            if (l < m - 1)
                varphi = ones(n, 3);
                varphi(:, 1) = calc_EuroPut(exp(yk(ii)) * samplePath0(j, :, l)', delta_t, r, sigma, K) - c1;
                varphi(:, 2) = calc_EuroPut(exp(yk(ii)) * samplePath0(j, :, l)', delta_t * (m - l), r, sigma, K) - c2;
                dd = 3;
            else
                varphi = ones(n, 2);
                varphi(:, 1) = calc_EuroPut(exp(yk(ii)) * samplePath0(j, :, l)', delta_t, r, sigma, K) - c1;
                dd = 2;
            end
            Psi_n = (varphi') * varphi / n;
            %yy = sampleY(:, j) .* coeff / rr;
            yy = sampleY(:, j);
            Phi_n = (yy') * varphi / n;
            Upsilon_n = (varphi') * (coeff .* signs) / n;
            ml = Phi_n / Psi_n;
            mr = Psi_n \ Upsilon_n;
            
            mat1 = repmat(reshape(varphi, n, dd, 1), [1, 1, dd]);
            mat2 = repmat(reshape(varphi, n, 1, dd), [1, dd, 1]);
            mat1 = mat1 .* mat2;
            mat1 = permute(mat1, [2 3 1]);
            mat1 = reshape(mat1, dd, dd * n);
            mat1 = ml * mat1;
            mat1 = reshape(mat1, dd, n);
            mat1 = (mat1') * mr;
            v2 = v2 + mat1;
            mat1 = repmat(yy, [1, dd]);
            mat1 = mat1 .* varphi;
            mat1 = mat1 * mr;
            v2 = v2 - mat1;
        end
    end
end



UB_sigma2 = sum((v2 .* v2), 'all') / n;
UB_sigma2 = UB_sigma2 - (sum(z .* signs, 'all') / n)^2;
UB_sigma2 = UB_sigma2 / ((1 - alpha)^2);

y = [LB_estimate, LB_sigma2, UB_estimate, UB_sigma2, LS_estimate];