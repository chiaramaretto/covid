function dx = SIRVaccineRecoveryFunction(t, x, beta, gamma, M, h,c,m, N, dist, C_max, omega, T_osp, T_icu, restrictions, T, ritardo)
% Versione con 
% 1) vaccini
% 2) restrizioni opzionali
% 3) perdita di immunità
% 4) svuotamento ospedale
    n = size(M, 1); %dovrebbe essere 9
    dx = zeros(6*n, 1);
    
    % Suddivisione del vettore di stato
    S = x(1:n);
    I = x(n+1:2*n);
    R = x(2*n+1:3*n);
    H = x(3*n+1:4*n);
    C = x(4*n+1:5*n);
    D = x(5*n+1:6*n);
  
    if ritardo
        if t > 300
            mu = T * (omega .* S); 
            mu = max(mu, 0);
        else
            mu = zeros(n, 1); 
        end
    else
        mu = T * (omega .* S); 
        mu = max(mu, 0);
    end
    

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
    dx(1:n) = -beta * (S ./ N) .* (M_change * I) + 1/200 * (R-H-D) - mu ; %modifica
    
    % dI_i/dt
    dx(n+1:2*n) = beta * (S ./ N) .* M_change * I - gamma * I;
    
    % dR_i/dt
    dx(2*n+1:3*n) = gamma * I - 1/200 * (R-H-D) + mu; % modifica

    % dH_i/dt
    dx(3*n+1:4*n) = gamma * h .* I - 1/T_osp * H;

    % dC_i/dt
    dx(4*n+1:5*n) = gamma * h .* c .* I - 1/T_icu * C;

    % dM_i/dt
    dx(5*n+1:6*n) = gamma * m .* I + 1.5/T_osp * m .* H + 2/T_icu * m .* C;
end