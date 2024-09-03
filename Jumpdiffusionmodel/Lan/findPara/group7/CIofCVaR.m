function y = CIofCVaR(D, S0, K, mu, r, T, sigma, J_lambda, J_m, J_theta,tau, l2,alpha, alphao, alphai,alphas,alphahi,alphalo,p,Budget,c,n0,m0,delta,lmin,lmax)
% CIofCVaR output CI = [L, U]  
CorMatrix = 0.3*ones(D,D)+0.7*diag(ones(1,D));
temp = chol(CorMatrix);
CorL = temp'; 
clear temp


V0 = calc_JumpEuroPut(S0, T, r, sigma, K, J_lambda, J_m, J_theta);
J_r1 = mu - (sigma^2)/2 - (exp(J_m + (J_theta^2)/2) - 1) * J_lambda;
SampleXn0 = S0*exp(J_r1*tau*ones(D,n0) + sigma*sqrt(tau)*CorL*normrnd(0,1,D,n0));
J_N01 = normrnd(0, 1, D, n0);
J_J = random('Poisson', tau * J_lambda, D, n0);
SampleXn0 = SampleXn0 .* exp(J_N01 * J_theta .* sqrt(J_J) + J_m * J_J);
SampleS_N = NestedEsti(sigma,J_lambda, J_m, J_theta,SampleXn0,m0,D,T,tau,l2,r); 
SampleYm0 = max(K - SampleS_N,0);
SampleYm0 = SampleYm0 * exp(-r*(T-tau));
clear SampleS_N;
SampleYm0 = SampleYm0 - V0;
SampleYm0 = squeeze(sum(SampleYm0, 1));
Ymeanm0 = mean(SampleYm0(:,:),2);
Yvarm0 = zeros(1,n0);
for i = 1:n0
    Yvarm0(i) = var(SampleYm0(i,:));
end
[SortYmean, Index0] = sort(Ymeanm0); % [3 1 2] becomes [1 2 3] Index = [2 3 1]
clear SortYmean;
d = tinv(1-alphas/((n0-round(n0*p))*round(n0*p)), m0-1);

I = zeros(1,n0);
I(1:round(n0*p)) = 1:round(n0*p);

i = n0;
while (i>round(n0*p))
    b = 0;
    j = 1;
    while (j<i)
        diffij = SampleYm0(Index0(i),:) - SampleYm0(Index0(j),:);
        S2m0ij = var(diffij);
        if (Ymeanm0(Index0(i))>Ymeanm0(Index0(j))+d*(S2m0ij)^(1/2)/sqrt(m0))
            b = b + 1;
        end
        j = j + 1;
    end
    if (b<round(n0*p))
        I(i) = i;
    end
    i = i - 1;
end 
I = I(I~=0); 
clear Ymeanm0 Yvarm0;

Budget1 = Budget - n0*m0;
S2m0 = zeros(1,length(I));
for i = 1:length(I)
    S2m0(i) = var(SampleYm0(Index0(I(i)),:));
end
M = round(Budget1*S2m0/sum(S2m0));
maxM = max(M);
clear sampleYm0;

SampleYS_N1 = NestedEsti(sigma,J_lambda, J_m, J_theta,SampleXn0(:, Index0(I)), maxM,D,T,tau,l2,r);
SampleYM = max(K - SampleYS_N1,0);
SampleYM = SampleYM * exp(-r*(T-tau));
clear sampleYS_N1 SampleXn0;
SampleYM = SampleYM - V0;
SampleYM = squeeze(sum(SampleYM, 1));

for i = 1:length(I)
    SampleYM(i,(M(i)+1):maxM) = zeros(1, maxM-M(i));
end


YmeanM = zeros(1,length(I));
YvarM = zeros(1,length(I));
s = zeros(1,length(I)); 
for i = 1:length(I)
    YmeanM(i) = mean(SampleYM(i,1:M(i)));
    YvarM(i) = var(SampleYM(i,1:M(i)));
    s(i) = sqrt(YvarM(i)/M(i));
end

hatL = Inf;
tlo = zeros(1,n0+1);
Mlo = zeros(1,n0+1);
Mlo(floor(n0*p)) = min(M(1:floor(n0*p)));
s_underline = zeros(1,n0+1);
s_underline(floor(n0*p)) = max(s(1:floor(n0*p)));
B0 = zeros(1,lmax+1);
lmax_ = lmax;
if length(I) < lmax_
    lmax_ = length(I);
end
for l = floor(n0*p):lmax_
   tlo(l) = tinv(1-alphalo, Mlo(l)-1);
   B0(l) = s_underline(l)*delta(l); 
   res = Optimization(l,p,n0,c,YmeanM,0);
   hatL = min(hatL, res - tlo(l)*B0(l));
   if (l < lmax_)
       Mlo(l+1) = min(Mlo(l), M(l+1));
       s_underline(l+1) = max(s_underline(l), s(l+1));
   end
end    

SortYmean1 = sort(YmeanM);

hatU = -Inf;
s_bar = max(s);
M_hi = min(M);
t_hi = tinv(1-alphahi,M_hi-1);
BS = zeros(1,n0);
for l = lmin:round(n0*p)
   BS(l) = s_bar*delta(l);
   
   res = Optimization(l,p,n0,c,SortYmean1,1);
   hatU = max(hatU, res + t_hi*BS(l));
end
y = [hatL hatU];


