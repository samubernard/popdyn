function xx = bifurcation(map,x0,par,n)


burnin = 1000;
ns = burnin + n;
np = length(par);

xx = zeros(np,n);

for i=1:length(par)
    x = feval(map,x0,par(i),ns);
    xx(i,:)=x(burnin+1:ns);
end

