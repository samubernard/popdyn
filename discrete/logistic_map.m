function x=logistic_map(x0,mu,n)
% LOGISTIC_MAP discrete logistic map
%   x=logistic_map(x0,mu,n) first N iterates of the
%   logistic map x(t+1)=MU*x(t)*(1-x(t)) with x(0)=x0

x=zeros(1,n);
x(1)=x0;

for i = 2:n
    x(i)=logistic(x(i-1),mu,1);
end

