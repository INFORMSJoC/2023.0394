function y = NestedEsti(sigma,SampleX,inN,D,T,tau,r)

CorMatrix = 0.3*ones(D,D)+0.7*diag(ones(1,D));
temp = chol(CorMatrix);
CorL = temp'; %%-- A lower triangular matrix
clear temp
 
temp1 = size(SampleX);
n = temp1(2); 
clear temp1
Y = zeros(D, n, inN);
for j = 1:inN
    Y(:, :,j) = SampleX.*exp((r-sigma^2/2)*(T-tau)*ones(D,n) + sigma*sqrt(T-tau)*CorL*normrnd(0,1,D,n));
end
clear n inN
y = Y;% D-by-n-by-inN

