function dx = SIRAgeFunction(t, x, beta, gamma, M, N)
    n = size(M, 1); %dovrebbe essere 9
    dx = zeros(3*n, 1);
    
    % Suddivisione del vettore di stato
    S = x(1:n);
    I = x(n+1:2*n);
    R = x(2*n+1:3*n);
    
    % dS_i/dt
    dx(1:n) = -beta * (S ./ N) .* M * I;
    
    % dI_i/dt
    dx(n+1:2*n) = beta * (S ./ N) .* M * I - gamma * I;
    
    % dR_i/dt
    dx(2*n+1:3*n) = gamma * I;

    % per verificare se la condizione R_0=1 è adeguata per prevedere il
    % picco
    R_0 = max(eig(diag(S./N)*M));
    if abs(R_0 - gamma/beta)<10^-1
        beta/gamma*R_0
        t
    end
        
end