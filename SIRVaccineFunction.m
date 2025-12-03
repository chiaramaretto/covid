function dx = SIRVaccineFunction(t, x, beta, gamma, M, N, T, h,c,m, omega, restrictions, C_max, dist)
    n = size(M, 1); %dovrebbe essere 9
    dx = zeros(6*n, 1);
    
    % Suddivisione del vettore di stato
    S = x(1:n);
    I = x(n+1:2*n);
    R = x(2*n+1:3*n);
    H = x(3*n+1:4*n);
    C = x(4*n+1:5*n);
    D = x(5*n+1:6*n);


    % CHIEDERE A PREZIOSI: PUò ESSERE CHE IL PROBLEMA SIA SOLO LA SCALA?
    % Sembra che diminuendo mu vada
    mu = T * (omega .* S); % / sum(omega.* S);
    %mu = min(mu, S);
    mu = max(mu, 0);
    %mu = 0.02 * S;

    if restrictions == true
        C_sum = sum(C);
        if C_sum < C_max/dist
            M_change = M;
        else
            M_change = 1/(dist*C_sum/C_max)*M;
        end
    else
        M_change = M;
    end    

    % dS_i/dt
    dx(1:n) = -beta * (S ./ N) .* (M_change * I) - mu ;
    
    % dI_i/dt
    dx(n+1:2*n) = beta * (S ./ N) .* M_change * I - gamma * I;
    
    % dR_i/dt
    dx(2*n+1:3*n) = gamma * I + mu;

    % dH_i/dt
    dx(3*n+1:4*n) = gamma * h .* I;

    % dC_i/dt
    dx(4*n+1:5*n) = gamma * h .* c .* I;

    % dM_i/dt
    dx(5*n+1:6*n) = gamma * m .* I;
end