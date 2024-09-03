function y = construct_CI(S0, D, n, m, l1, T, r, mu, sigma, K, H, alpha)
CorMatrix = 0.3 * ones(D, D) + 0.7 * diag(ones(1, D));
temp = chol(CorMatrix);
CorL = temp';
V0 = calc_V0(S0, T, r, sigma, K, H);
delta_t = T / m;
tau = l1 * delta_t;
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
sampleMin = zeros(D, n, m);
UnifRand = rand(D, n, m);
sampleMin(:, :, 1) = exp((log(S0 * samplePath(:, :, 1)) - sqrt(power(log(samplePath(:, :, 1) / S0), 2) - 2 * sigma^2 * delta_t * log(UnifRand(:, :, 1))))/2);
for l = 2:m
    sampleMin(:, :, l) = exp((log(samplePath(:, :, l-1) .* samplePath(:, :, l)) - sqrt(power(log(samplePath(:, :, l) ./ samplePath(:, :, l-1)), 2) - 2 * sigma^2 * delta_t * log(UnifRand(:, :, l))))/2);
end
sampleX = zeros(2 * D, n);
sampleX(1:D, :) = samplePath(:, :, l1);
sampleX(D+1:D*2, :) = min(sampleMin(:, :, [1 : l1]), [], 3);

sampleY = max(samplePath(:, :, m) - K(1), 0) + max(samplePath(:, :, m) - K(2), 0) + max(samplePath(:, :, m) - K(3), 0);
sampleY = sampleY * exp(-r * (T - tau));
pathMin = min(sampleMin(:, :, :), [], 3);
sampleY = sampleY .* (pathMin >= H);
sampleY = ones(D, n) * V0 - sampleY;
sampleY = sampleY';
clear pathMin;

for l = l1:(m - 1)
    m0 = samplePath(:, :, l)';
    m0 = bsxfun(@rdivide, bsxfun(@minus, m0, mean(m0)), std(m0)); %(n)*(2D)
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
        minPre = (min(sampleMin(j, :, [1 : l]), [], 3))';
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
        ybar = varphi * beta .* (minPre >= H);
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
clear mat1 mat2 samplePath
UB_sigma2 = sum((v2 .* v2), 'all') / n;
UB_sigma2 = UB_sigma2 - (sum(z .* signs, 'all') / n)^2;
UB_sigma2 = UB_sigma2 / ((1 - alpha)^2);

sampleY = sum(sampleY, 2);
m0 = sampleX(:, :)';
m0 = bsxfun(@rdivide, bsxfun(@minus, m0, mean(m0)), std(m0));
varphi = ones(n, 8 * D + 1);
varphi(:, 1 : D*2) = exp(-m0 / 2.0);
varphi(:, D*2+1 : D*4) = varphi(:, 1 : D*2) .* (1 - m0);
varphi(:, D*4+1 : D*6) = varphi(:, D*2+1 : D*4) .* (1 - 2 * m0 + (m0 .^ 2) / 2);
varphi(:, D*6+1 : D*8) = varphi(:, D*4+1 : D*6) .* (1 - 3 * m0 + 3 * (m0 .^ 2) / 2 - (m0 .^ 3) / 6);
varphi(:, D*8+1) = 1;
clear m0 sampleX;
beta = varphi \ sampleY; 
ybar = varphi * (beta);
aveY = sum(sampleY, 'all') / n;
sortedY = [ybar, sampleY];
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