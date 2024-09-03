function y = NestedEsti(sigma,J_lambda, J_m, J_theta, SampleX,inN,D,T,tau,l2,r)

CorMatrix = 0.3*ones(D,D)+0.7*diag(ones(1,D));
temp = chol(CorMatrix);
CorL = temp'; %%-- A lower triangular matrix
clear temp
 
temp1 = size(SampleX);
n = temp1(2); 
clear temp1
Y = zeros(D, n, inN);
delta_t = (T-tau)/l2;
for j=1:inN
    Y(:,:,j) = SampleX;
end

J_r2 = r - (sigma^2)/2 - (exp(J_m + (J_theta^2)/2) - 1) * J_lambda;
for xx = 1:l2
    for j = 1:inN
        J_N01 = normrnd(0, 1, D, n);
        J_J = (rand([D, n]) >= exp(-J_lambda * delta_t));
        Y(:, :,j) = Y(:, :,j).*exp(J_r2*delta_t*ones(D,n) + sigma*sqrt(delta_t)*CorL*normrnd(0,1,D,n));
        Y(:,:,j) = Y(:,:,j) .* exp(J_N01 * J_theta .* J_J + J_m * J_J);
    end
end
clear n inN
y = Y;% D-by-n-by-inN

