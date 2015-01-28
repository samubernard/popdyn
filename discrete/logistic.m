function y = logistic(x,mu,p)
% LOGISTIC is the logistic function mu*x*(1-x^2)
%   y = logistic(x,mu,p) is the P-th iterate of the function
%   MU*X*(1-X)

if nargin == 2
    p = 1
end

y = x;
for i = 1:p
    y = f(y);
end


    function y = f(x)

        y = mu*x.*(1-x);

    end

end
