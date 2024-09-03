function y = Optimization(l,p,n0,c,YmeanM,direction)
x= sym('x',[1,2]);
eqn1 = 0;
temp = 0;
for i = 1:l
    eqn1 = eqn1 + (YmeanM(i)+x(1))/(1-(YmeanM(i)+x(1))*x(2));
    temp = temp + log(1-(YmeanM(i)+x(1))*x(2));
end 
eqn2 = n0*log(n0) - log(c)+(n0-l)*log((1-p)/(n0-l)) + l*log(p/l)-temp;
f= [eqn1,eqn2];
func = matlabFunction(f,'Vars',{[x(1), x(2)]});
options = optimoptions('fsolve','Display','off');
if direction == 0
    solution = fsolve(func,[-YmeanM(l),1],options);
else
    solution = fsolve(func,[-YmeanM(1),1],options);
end
y = real(solution(1));

