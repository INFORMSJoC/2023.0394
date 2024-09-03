function y = MAX(l,p,c,n0)
Max = zeros(1,l-1);

for m = 1:(l-1) 
    h = (l-m)*log(l-m)+ log(c)-n0*log(n0)-(n0-l)*log(1-p)+(n0-l)*log(n0-l);
    fmax = log(c)-n0*log(n0)-(n0-l)*log(1-p)+(n0-l)*log(n0-l)-l*log(p) +l*log(l);
    if (fmax==0)
        a = 1/l;
        Max(m) = m*a^2 + (l-m)*((1-m*a)/(l-m))^2;
    elseif (fmax>0)
        syms a;
        eq1 = m*log(a)+(l-m)*log(1-m*a)-h; 
        a = solve(eq1,a,'Real',true); 
        Max(m) = max( m*a.^2 + (l-m)*((1-m*a)/(l-m)).^2 );
    end
end

y = max(Max);

