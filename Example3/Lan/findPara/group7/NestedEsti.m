function y = NestedEsti(sigma,SampleX,inN,D,delta_t,r,l2)

CorMatrix = 0.3*ones(D,D)+0.7*diag(ones(1,D));
temp = chol(CorMatrix);
CorL = temp'; %%-- A lower triangular matrix
clear temp
 
temp1 = size(SampleX);
n = temp1(2); 
clear temp1
samplePath0 = repmat(reshape(SampleX, 2*D,n,1),[1,1,inN]);
samplePath1 = zeros(D*2, n, inN);
for l = 1:l2
    for j = 1:inN
        samplePath1(1:D, :, j) = exp((r - sigma^2 / 2) * delta_t ...
            + sigma * CorL * sqrt(delta_t) * normrnd(0, 1, D, n)) .* samplePath0(1:D, :, j);
        samplePath1(D+1:D*2, :, j) = exp((log(samplePath0(1:D, :, j) .* samplePath1(1:D, :, j)) - sqrt(power(log(samplePath1(1:D, :, j) ./ samplePath0(1:D, :, j)), 2) - 2 * sigma^2 * delta_t * log(rand(D, n))))/2);
        samplePath1(D+1:D*2, :, j) = min(samplePath0(D+1:D*2, :, j), samplePath1(1:D, :, j));
    end
    samplePath0 = samplePath1;
end
y = samplePath1;

