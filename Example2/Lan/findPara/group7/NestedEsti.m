function y = NestedEsti(sigma,SampleX,inN,D,delta_t,r,l2, m)

CorMatrix = 0.3*ones(D,D)+0.7*diag(ones(1,D));
temp = chol(CorMatrix);
CorL = temp'; %%-- A lower triangular matrix
clear temp
 
temp1 = size(SampleX);
n = temp1(2); 
clear temp1
Y0 = repmat(reshape(SampleX, 2*D,n,1),[1,1,inN]);
Y1 = zeros(D*2, n, inN);
for l = 1:l2
    for j = 1:inN
        Y1(1:D, :, j) = Y0(1:D, :, j) + (r-sigma^2/2)*delta_t ...
            + sigma*sqrt(delta_t)*CorL*normrnd(0,1,D,n);
        Y1(D+1:2*D, :, j) = Y0(D+1:2*D, :, j) + Y1(1:D, :, j);
    end
    Y0 = Y1;
end
y = exp(Y1(D+1:2*D, :, :) / m);

