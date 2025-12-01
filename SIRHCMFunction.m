function dx = SIRHCMFunction(t, x, beta, gamma, M, h,c,m, N)
    n = size(M, 1); %dovrebbe essere 9
    dx = zeros(6*n, 1);
    
    % Suddivisione del vettore di stato
    S = x(1:n);
    I = x(n+1:2*n);
    R = x(2*n+1:3*n);
    H = x(3*n+1:4*n);
    C = x(4*n+1:5*n);
    D = x(5*n+1:6*n);

    
    % dS_i/dt
    dx(1:n) = -beta * (S ./ N) .* (M * I);
    
    % dI_i/dt
    dx(n+1:2*n) = beta * (S ./ N) .* M * I - gamma * I;
    
    % dR_i/dt
    dx(2*n+1:3*n) = gamma * I;

    % dH_i/dt
    dx(3*n+1:4*n) = gamma * h .* I;

    % dC_i/dt
    dx(4*n+1:5*n) = gamma * h .* c .* I;

    % dM_i/dt
    dx(5*n+1:6*n) = gamma * m .* I;
end