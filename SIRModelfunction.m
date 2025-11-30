function dx = SIRModelfunction(t, x, beta, gamma)
    dx= zeros(3,1);
    dx(1) = -beta*x(1)*x(2);
    dx(2) = beta*x(1)*x(2)-gamma*x(2);
    dx(3) = gamma*x(2);
end