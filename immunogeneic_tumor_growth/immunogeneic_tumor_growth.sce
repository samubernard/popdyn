// IMMUNOGENEIC_TUMOR_GROWTH model for tumor immune interaction
//
//   REFERENCE: Kuznetsov et al. (1994) Bull Math Biol 56, 295-321
//
//   University Lyon 1 - 2014-2015

function dXdt = oderhs(t,X)
    // F rhs of the ODE system dX/dt = f(X)
    //
    //   X = (x,y)
    //
    //   dx/dt = sigma + rho*x*y/(eta + y) - mu*x*y - delta*x
    //   dy/dt = alpha*y*(1-beta*y) - x*y


    // Non-dimensional model 
    // System parameters
    sigma = 0.1181;
    rho = 1.131;
    eta = 20.19;
    mu = 0.00311;
    delta = 0.3743;
    alpha = 1.636;
    beta_par = 1e-2; 

    x = X(1);
    y = X(2);

    dXdt = [sigma + rho*x.*y./(eta + y) - mu*x.*y - delta*x;
            alpha*y.*(1-beta_par*y) - x.*y];

endfunction

function x = nullx(y)
    // NULLX nullcline for x (dxdt = 0 when x = nullx(y))

    // Non-dimensional model 
    // System parameters
    sigma = 0.1181;
    rho = 1.131;
    eta = 20.19;
    mu = 0.00311;
    delta = 0.3743;
    alpha = 1.636;
    beta_par = 1e-2; 

    x = sigma./(delta+mu*y-rho*y./(eta+y));

endfunction

function x = nully(y)
    // NULLX nullcline for y (dydt = 0 when x = nully(y))

    // Non-dimensional model 
    // System parameters
    sigma = 0.1181;
    rho = 1.131;
    eta = 20.19;
    mu = 0.00311;
    delta = 0.3743;
    alpha = 1.636;
    beta_par = 1e-2; 

    x = alpha*(1-beta_par*y);

endfunction

// Simulations
// Compute a single trajectory with defined initial condition in interval [0,tfinal]
ic = [0.7616; 268.7980];
t0 = 0;
tspan = 0:100;
X = ode(ic,t0,t,oderhs);

//// Plot the trajectory in the phase space
figure(1); clf;
plot(X(1,:),X(2,:))
// plot the positive x-nullcline
yy = logspace(-1,3,300);
yyp = yy(nullx(yy)>0)
plot(nullx(yyp),yyp,'k')
// plot the y-nullcline
plot(nully(yy),yy,'r')
mtlb_axis([1e-1 5 1e-1 1e3])
xlabel('x')
ylabel('y')








