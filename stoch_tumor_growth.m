function [tt,x,y] = stoch_tumor_growth(x0,y0)

%% Non-dimensional model 
%% System parameters
sigma = 0.1181;
rho = 1.131;
eta = 20.19;
mu = 0.00311;
delta = 0.3743;
alpha = 1.636;
beta = 2.0e-3;

k=1;
t=0;
x(k)=x0;
y(k)=y0;
tt(k)=t;

tfinal = 100;


while t<tfinal
    r = [ sigma;
          rho*x(k)*y(k)/(eta+y(k));
          mu*x(k)*y(k);
          delta*x(k);
          alpha*y(k);
          alpha*beta*y(k)*y(k);
          x(k)*y(k)];
    sum_r = sum(r);
    dt = exprnd(1/sum_r); % exprnd takes the mean waiting time as argument, not the rate
    u = rand(1);
    c = cumsum(r/sum_r);
    idx = min(find(u < c));
    switch idx
        case 1 % birth of x
            x(k+1)=x(k)+1;
            y(k+1)=y(k);
        case 2 % birth of x
            x(k+1)=x(k)+1;
            y(k+1)=y(k);
        case 3 % death of x
            x(k+1)=x(k)-1;
            y(k+1)=y(k);
        case 4 % death of x
            x(k+1)=x(k)-1;
            y(k+1)=y(k);
        case 5 % birth of y
            y(k+1)=y(k)+1;
            x(k+1)=x(k);
        case 6 % death of y
            y(k+1)=y(k)-1;
            x(k+1)=x(k);
        case 7 % death of y
            y(k+1)=y(k)-1;
            x(k+1)=x(k);
    end
    t = t+dt;
    tt(k+1)=t;
    k=k+1;
    %fprintf('t=%f, x=%f, y=%f\n',t,x(k),y(k));

    if y(k)==0
        fprintf('tumor eradicated\n');
        break        
    end
end

