function y = construct_CI(S0, D, n, T, m, l1, r, mu, sigma, K, alpha)

CorMatrix = 0.3 * ones(D, D) + 0.7 * diag(ones(1, D));
temp = chol(CorMatrix);
CorL = temp';
V0 = calc_V0(S0, T, r, sigma, K);
delta_t = T / m;
tau = delta_t * l1;
N01 = normrnd(0, 1, D, n, m);
samplePath = zeros(D, n, m);
samplePath(:, :, 1) = exp((mu - sigma^2 / 2) * delta_t * ones(D, n) ...
    + sigma * CorL * sqrt(delta_t) * N01(:, :, 1)) * S0;

for l = 2:l1
    samplePath(:, :, l) = exp((mu - sigma^2 / 2) * delta_t * ones(D, n) ...
    + sigma * CorL * sqrt(delta_t) * N01(:, :, l)) .* samplePath(:, :, l - 1);
end
for l = (l1 + 1):m
    samplePath(:, :, l) = exp((r - sigma^2 / 2) * delta_t * ones(D, n) ...
    + sigma * CorL * sqrt(delta_t) * N01(:, :, l)) .* samplePath(:, :, l - 1);
end
tmp = samplePath(:, :, m);
tmp = max(tmp - K(1), 0.0) + max(tmp - K(2), 0.0) + max(tmp - K(3), 0.0);
sampleY = V0 * ones(D, n) - tmp / exp((T - tau) * r);
sampleY = sampleY'; % (n)*(D)
clear tmp;

for l = l1:(m - 1)
    m0 = (samplePath(:, :, l))';
    m0 = bsxfun(@rdivide, bsxfun(@minus, m0, mean(m0)), std(m0)); %(n)*(D)
    coBM = sqrt(delta_t) * CorL * N01(:, :, l + 1);
    samplePath(:, :, l) = m0';
    N01(:, :, l + 1) = coBM;
end

delta_martingale = zeros(n, 1);

varphi = zeros(n, 5);
Psi_n = zeros(5, 5, m, m);
Phi_n = zeros(1, 5, m, m);
for l = l1:(m - 1)
    for j = 1:D
        m1 = samplePath(j, :, l)';
        varphi(:, 1) = exp(-m1 / 2.0);
        varphi(:, 2) = varphi(:, 1) .* (1 - m1);
        varphi(:, 3) = varphi(:, 2) .* (1 - 2 * m1 + (m1 .^ 2) / 2);
        varphi(:, 4) = varphi(:, 3) .* (1 - 3 * m1 + 3 * (m1 .^ 2) / 2 - (m1 .^ 3) / 6);
        varphi(:, 5) = 1;
        coBM_j = N01(j, :, l + 1)';
        beta = varphi \ sampleY(:, j);
        y_residual = sampleY(:, j) - varphi * beta;
        beta = varphi \ (y_residual .* coBM_j / delta_t);
        ybar = varphi * beta;
        delta_martingale = delta_martingale + ybar .* coBM_j;
        Psi_n(:, :, l, j) = (varphi') * varphi / n;
        yy = sampleY(:, j) .* coBM_j / delta_t;
        Phi_n(:, :, l, j) = (yy') * varphi / n;
    end
end
z = sum(sampleY, 2) - delta_martingale;
z = sort(z, 'descend');
sz = ceil((1 - alpha) * n);
UB_u = z(sz);
UB_estimate = sum(z(1 : sz), 'all') / sz;

z = sum(sampleY, 2) - delta_martingale - UB_u;
signs = (z >= 0);
Upsilon_n = zeros(5, 1, m, m);
for l = l1:(m - 1)
    for j = 1:D
        m1 = samplePath(j, :, l)';
        varphi(:, 1) = exp(-m1 / 2.0);
        varphi(:, 2) = varphi(:, 1) .* (1 - m1);
        varphi(:, 3) = varphi(:, 2) .* (1 - 2 * m1 + (m1 .^ 2) / 2);
        varphi(:, 4) = varphi(:, 3) .* (1 - 3 * m1 + 3 * (m1 .^ 2) / 2 - (m1 .^ 3) / 6);
        varphi(:, 5) = 1;
        coBM_j = N01(j, :, l + 1)';
        Upsilon_n(:, :, l, j) = (varphi') * (coBM_j .* signs) / n;
    end
end
ml = zeros(1, 5, m, m);
mr = zeros(5, 1, m, m);
for l = l1:(m-1)
    for j = 1:D
        ml(:, :, l, j) = Phi_n(:, :, l, j) / Psi_n(:, :, l, j);
        mr(:, :, l, j) = Psi_n(:, :, l, j) \ Upsilon_n(:, :, l, j);
    end
end
v2 = z .* signs;
for l = l1:(m - 1)
    for j = 1:D
        m1 = samplePath(j, :, l)';
        varphi(:, 1) = exp(-m1 / 2.0);
        varphi(:, 2) = varphi(:, 1) .* (1 - m1);
        varphi(:, 3) = varphi(:, 2) .* (1 - 2 * m1 + (m1 .^ 2) / 2);
        varphi(:, 4) = varphi(:, 3) .* (1 - 3 * m1 + 3 * (m1 .^ 2) / 2 - (m1 .^ 3) / 6);
        varphi(:, 5) = 1;
        coBM_j = N01(j, :, l + 1)';
        yy = sampleY(:, j) .* coBM_j / delta_t;
        
        mat1 = repmat(reshape(varphi, n, 5, 1), [1, 1, 5]);
        mat2 = repmat(reshape(varphi, n, 1, 5), [1, 5, 1]);
        mat1 = mat1 .* mat2; %(n)*(5)*(5)
        mat1 = permute(mat1, [2 3 1]);
        mat1 = reshape(mat1, 5, 5 * n);
        mat1 = ml(:, :, l, j) * mat1;
        mat1 = reshape(mat1, 5, n);
        mat1 = (mat1') * mr(:, :, l, j);
        v2 = v2 + mat1;
        mat1 = repmat(yy, [1, 5]);
        mat1 = mat1 .* varphi;
        mat1 = mat1 * mr(:, :, l, j);
        v2 = v2 - mat1;
    end
end
clear mat1 mat2
UB_sigma2 = sum((v2 .* v2), 'all') / n;
UB_sigma2 = UB_sigma2 - (sum(z .* signs, 'all') / n)^2;
UB_sigma2 = UB_sigma2 / ((1 - alpha)^2);

sampleY = sum(sampleY, 2);
m0 = samplePath(:, :, l1)';
clear samplePath;
varphi = ones(n, 4 * D + 1);
varphi(:, 1 : D) = exp(-m0 / 2.0);
varphi(:, D+1 : D*2) = varphi(:, 1 : D) .* (1 - m0);
varphi(:, D*2+1 : D*3) = varphi(:, D+1 : D*2) .* (1 - 2 * m0 + (m0 .^ 2) / 2);
varphi(:, D*3+1 : D*4) = varphi(:, D*2+1 : D*3) .* (1 - 3 * m0 + 3 * (m0 .^ 2) / 2 - (m0 .^ 3) / 6);
varphi(:, D*4+1) = 1;
%range of elements in varphi: about 1-80
clear m0;
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

Psi_n = (varphi') * varphi / n;
Phi_n = (sampleY') * varphi / n;
epsilon = std(ybar) * 1.06 * n^(-1/5);
Upsilon_n = (sampleY - LB_u) .* (1 - cos(ybar - LB_u)) .* (abs(ybar - LB_u) <= 2 * pi * epsilon) / (4 * pi * epsilon);
Upsilon_n = ((Upsilon_n') * varphi)' / n;
v2 = 0;
m1 = Phi_n / Psi_n;
m2 = Psi_n \ Upsilon_n;

for k = 1:n
    Psi_k = varphi(k,:);
    Psi_k = (Psi_k') * (Psi_k);
    Phi_k = sampleY(k, 1) * varphi(k, :);
    v = -m1 * Psi_k * m2;
    v = v + Phi_k * m2;
    v = v + (sampleY(k, 1) - LB_u) * (ybar(k, 1) >= LB_u);
    v2 = v2 + v * v;
end
LB_sigma2 = v2 / n - (sum((sampleY - LB_u) .* (ybar >= LB_u), 'all') / n) ^ 2;
LB_sigma2 = LB_sigma2 / ((1 - alpha)^2);


%LB_estimate
%LB_sigma2
%UB_estimate
%UB_sigma2
y = [LB_estimate, LB_sigma2, UB_estimate, UB_sigma2, LS_estimate];
