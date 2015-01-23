function out = immunogeneic_tumor_growth()
% IMMUNOGENEIC_TUMOR_GROWTH model for tumor immune interaction
%   out = immunogeneic_tumor_growth()
%
%   REFERENCE: Kuznetsov et al. (1994) Bull Math Biol 56, 295-321
%
%   Université Lyon 1 - 2014-2015

%% Non-dimensional model 
%% System parameters
sigma = 0.1181;
rho = 1.131;
eta = 20.19;
mu = 0.00311;
delta = 0.3743;
alpha = 1.636;
beta = 1e-2; %2.0e-3;

%% Simulations
%% Compute a single trajectory with defined initial condition in interval [0,tfinal]
ic = [0.7616,268.7980];
tspan = [0,100];
options = odeset('RelTol',1e-3,'AbsTol',1e-6);
sol = ode45(@f,tspan,ic,options);

%% Plot the trajectory in the phase space
t = tspan(1):0.1:tspan(2);
X = deval(sol,t);
figure(1); clf;
loglog(X(1,:),X(2,:))
hold on;
% plot the positive x-nullcline
yy = logspace(-1,3,300);
yyp = yy(nullx(yy)>0)
loglog(nullx(yyp),yyp,'k')
% plot the y-nullcline
loglog(nully(yy),yy,'r')
axis([1e-1 5 1e-1 1e3])
xlabel('x')
ylabel('y')

%% Steady states
C0 = eta*(sigma/alpha - delta);
C1 = sigma/alpha + rho - mu*eta - delta + delta*eta*beta;
C2 = -mu + (mu*eta + delta - rho)*beta;
C3 = mu*beta;

C = [C3 C2 C1 C0];

stst_y = roots(C);
stst_x = sigma./(delta+mu*stst_y-rho*stst_y./(eta+stst_y));

stst_x =[stst_x; sigma/delta];
stst_y =[stst_y; 0];

stst = [stst_x, stst_y];
stst = stst(stst_x >= 0 & stst_y >= 0,:);

[nbr_stst,~] = size(stst);

%% Stability analysis
for i = 1:nbr_stst
    xx = stst(i,1);
    yy = stst(i,2);
    J = [rho*yy/(eta+yy) - mu*yy - delta, ((eta+yy)*rho*xx-rho*xx*yy)/(eta+yy)^2-mu*xx;
    -yy, alpha*(1-2*beta*yy)-xx];
    [V,D]=eig(J);
    fprintf('STST %d, (x,y) = (%f,%f)\n',i,xx,yy);
    J
    D
    V 
end



%% dynamical system definition
    
    function dXdt = f(~,X)
    % F rhs of the ODE system dX/dt = f(X)
    %
    %   X = (x,y)
    %
    %   dx/dt = sigma + rho*x*y/(eta + y) - mu*x*y - delta*x
    %   dy/dt = alpha*y*(1-beta*y) - x*y
    %
    %   Note: F is a nested function. It has access to the local variables 
    %   of its parent function immunogeneic_tumor_growth. Parameters sigma, rho, etc
    %   do not need to be passed to the nested function.


    x = X(1);
    y = X(2);

    dXdt = [sigma + rho*x.*y./(eta + y) - mu*x.*y - delta*x;
            alpha*y.*(1-beta*y) - x.*y];

    end

    function x = nullx(y)
    % NULLX nullcline for x (dxdt = 0 when x = nullx(y))

    x = sigma./(delta+mu*y-rho*y./(eta+y));

    end

    function x = nully(y)
    % NULLX nullcline for y (dydt = 0 when x = nully(y))

    x = alpha*(1-beta*y);

    end


end
