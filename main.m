clear all
close all
clc

%% Dati italiani
% https://esploradati.istat.it/databrowser/#/it/censpop/dashboards
% Sul sito istat ci sono più fasce d'età, con ampiezza 5 anni (fino a 99,
% poi 100+)
P = [2031476; 2409672; 2745119; 2921261; 2950359; 3012211; 3213337; ...
    3357101; 3680853; 4380844; 4770527; 4848988; 4292746; ...
    3672071; 3254150; 2887393; 2223335; 1481572; 654100; 162904; 21211];
N = sum(P); % controllo anche che sia coerente col totale istat OK
T_osp = 14.8; % tempo medio di degenza in ospedale
gamma = 1/14; % teniamo lo stesso del paper
T_icu = 9; 

%% ======================================================
%            Simulazione del modello SIR
% Model parameters
R_0 = 5; %5.7
M = [19.2, 4.8, 3.0, 7.1, 3.7, 3.1, 2.3, 1.4, 1.4;
     4.8, 42.4, 6.4, 5.4, 7.5, 5.0, 1.8, 1.7, 1.7;
     3.0, 6.4, 20.7, 9.2, 7.1, 6.3, 2.0, 0.9, 0.9;
     7.1, 5.4, 9.2, 16.9, 10.1, 6.8, 3.4, 1.5, 1.5;
     3.7, 7.5, 7.1, 10.1, 13.1, 7.4, 2.6, 2.1, 2.1;
     3.1, 5.0, 6.3, 6.8, 7.4, 10.4, 3.5, 1.8, 1.8;
     2.3, 1.8, 2.0, 3.4, 2.6, 3.5, 7.5, 3.2, 3.2;
     1.4, 1.7, 0.9, 1.5, 2.1, 1.8, 3.2, 7.2, 7.2;
     1.4, 1.7, 0.9, 1.5, 2.1, 1.8, 3.2, 7.2, 7.2];

n = size(M, 1);

% Popolazione italiana (per fasce 0-9, ..., 80+)
P = [910147; 963765; 1049452; 1141263; 959404; 951386; 904071; 561572; 265250];
N = sum(P);

% Calcolo proporzioni
p = P / N;
lambda = max(abs(eigs(M*diag(p))));
gamma = 1/14;
beta = R_0*gamma/lambda;

T = 365*2; 
S0 = P;
I0 = ones(n,1);
R0 = zeros(n,1);

age_groups = {'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
colors = lines(n);

h = [0.001; 0.003; 0.012; 0.032; 0.049; 0.102; 0.166; 0.243; 0.273];
c = [0.050; 0.050; 0.050; 0.050; 0.063; 0.122; 0.274; 0.432; 0.709];
m = [0.00002; 0.00006; 0.00030; 0.00080; 0.00150; 0.00600; 0.02200; 0.05100; 0.09300];

%% ======================================================
%                   SIR BASE con fasce d'età

clc

X0 = [S0; I0 ; R0];

[t,X] = ode45(@(t,x)SIRAgeFunction(t, x, beta, gamma, M, N),[0, T],X0);
S = X(:,1:n);   I = X(:, n+1:2*n);   R = X(:, 2*n+1:3*n);

age_groups = {'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
colors = lines(n);

% --- Grafici S, I, R per fasce ---
figure;
tiledlayout(3,1);

nexttile; hold on;
for k=1:n, plot(t,S(:,k)/P(k),'LineWidth',2,'Color',colors(k,:)); end
title('Susceptible by Age Group'); ylabel('S'); grid on;

nexttile; hold on;
for k=1:n, plot(t,I(:,k)/P(k),'LineWidth',2,'Color',colors(k,:)); end
title('Infected by Age Group'); ylabel('I'); grid on;

nexttile; hold on;
for k=1:n, plot(t,R(:,k)/P(k),'LineWidth',2,'Color',colors(k,:)); end
title('Removed by Age Group'); ylabel('R'); xlabel('Days'); grid on;

legend(age_groups,'Location','SouthOutside','NumColumns',5);

for k=1:9 % controllo 
    P(k) - (S(end,k) + I(end,k) + R(end,k)); 
    % giusto che venga -1 perché partiamo con S=P ma anche con Ii=1 per
    % ogni i
end

%% ======================================================
%               SIR con H,C,D sottoinsiemi di R
%               Parametri H, C, D (Ospedale, ICU, Morti)

clc

h = [0.001; 0.003; 0.012; 0.032; 0.049; 0.102; 0.166; 0.243; 0.273];
c = [0.050; 0.050; 0.050; 0.050; 0.063; 0.122; 0.274; 0.432; 0.709];
m = [0.00002; 0.00006; 0.00030; 0.00080; 0.00150; 0.00600; 0.02200; 0.05100; 0.09300];

H0 = zeros(n,1);
C0 = zeros(n,1);
M0 = zeros(n,1);
X0 = [S0; I0 ; R0; H0; C0; M0];

[t,X] = ode45(@(t,x)SIRHCMFunction(t, x, beta, gamma, M,h,c,m, N),[0, T],X0);

S = X(:,1:n); 
I = X(:, n+1:2*n);
R = X(:, 2*n+1:3*n);
H = X(:, 3*n+1:4*n);
C = X(:, 4*n+1:5*n);
D = X(:, 5*n+1:6*n);

% --- Grafici H, C, M per fasce ---
figure;
tiledlayout(3,1);

nexttile; hold on;
for k=1:n, plot(t,H(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Hospitalized H(t)'); ylabel('H'); grid on;

nexttile; hold on;
for k=1:n, plot(t,C(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Critical Care C(t)'); ylabel('C'); grid on;

nexttile; hold on;
for k=1:n, plot(t,D(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Deaths M(t)'); ylabel('M'); xlabel('Days'); grid on;

legend(age_groups,'Location','SouthOutside','NumColumns',5);

figure;
tiledlayout(3,1);

nexttile; hold on;
for k=1:n, plot(t,S(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Susceptible S(t)'); ylabel('S'); grid on;

nexttile; hold on;
for k=1:n, plot(t,I(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Infected I(t)'); ylabel('I'); grid on;

nexttile; hold on;
for k=1:n, plot(t,R(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Removed R(t)'); ylabel('R'); xlabel('Days'); grid on;

legend(age_groups,'Location','SouthOutside','NumColumns',5);

for k=1:9
    P(k) - (S(end,k) + I(end,k) + R(end,k)); % ok
end

%% ======================================================
%       Modello con restrizioni su M in base alle ICU

clc

% NOTA: METTI DEI VALORI SENSATI QUI (E FAI VARIARE dist)
dist = 0.5; % parametro che regola il livello di restrizioni applicate (più è alto, maggiori le restrizioni)
C_max = 34.7 * N /100000;

[t,X] = ode45(@(t,x)SIRHCMFunctionRestriction(t, x, beta, gamma, M,h,c,m, N, dist, C_max),[0, T],X0);

S = X(:,1:n); 
I = X(:, n+1:2*n);
R = X(:, 2*n+1:3*n);
H = X(:, 3*n+1:4*n);
C = X(:, 4*n+1:5*n);
D = X(:, 5*n+1:6*n);

% --- Grafici restrizioni ---
figure;
tiledlayout(3,1);

nexttile; hold on;
for k=1:n, plot(t,S(:,k),'LineWidth',2,'Color',colors(k,:)); end
title('Susceptible with Restrictions'); grid on;

nexttile; hold on;
for k=1:n, plot(t,I(:,k),'LineWidth',2,'Color',colors(k,:)); end
title('Infected with Restrictions'); grid on;

nexttile; hold on;
for k=1:n, plot(t,C(:,k),'LineWidth',2,'Color',colors(k,:)); end
title('ICU with Restrictions'); xlabel('Days'); grid on;

legend(age_groups,'Location','SouthOutside','NumColumns',5);

for k=1:9
    P(k) - (S(end,k) + I(end,k) + R(end,k)); %ok
end

%% ======================================================
%               Vaccinazione SENZA restrizioni
clc

X0 = [S0; I0 ; R0; H0; C0; M0];

vax_per_day = 1/720;

omega_C = [1; 1; 1; 1; 1; 1; 1; 1; 1]; % vaccinazioni uniformi per età
omega_M = [1; 1; 1; 1; 1; 2; 2; 4; 4]; % campagna di vaccinazioni con intensificazione per fasce fragili
omega_S = [1; 1; 1; 1; 2; 4; 8; 16; 16]; % campagna di vaccinazioni intensiva con intensificazione per fasce

omegas = {omega_C, omega_M, omega_S};

age_groups = {'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
colors = lines(n);
figure;
tiledlayout(3,3);  

for i = 1:3
    omega = omegas{i};

    [t,X] = ode45(@(t,x) SIRVaccineFunction(t, x, beta, gamma, M, N, vax_per_day, h,c,m, omega, false, [], []),[0, T],X0);

    S = X(:,1:n); 
    I = X(:, n+1:2*n);
    R = X(:, 2*n+1:3*n);
    H = X(:, 3*n+1:4*n);
    C = X(:, 4*n+1:5*n);
    D = X(:, 5*n+1:6*n);

    % versione dei grafici normalizzata
    % --- S(t) ---
    nexttile; hold on;
    for k=1:n, plot(t,S(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
    title(['Susceptible - \omega_' num2str(i)]); grid on;

    % --- I(t) ---
    nexttile; hold on;
    for k=1:n, plot(t,I(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
    title(['Infected - \omega_' num2str(i)]); grid on;

    % --- ICU C(t) ---
    nexttile; hold on;
    for k=1:n, plot(t,C(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
    title(['ICU - \omega_' num2str(i)]); grid on;

    norm(P-(S(end,:)+I(end,:)+R(end,:)))
    for k=1:9
        P(k) - (S(end,k) + I(end,k) + R(end,k))
        P(k) - R(end,k) - S(end, k)
    end
    
end
legend(age_groups,'Location','SouthOutside','NumColumns',5);




%% ======================================================
%               Vaccinazione CON restrizioni

vax_per_day = 1/720;

dist = 5; 
C_max = 34.7 * N /100000;

omega_C = [1; 1; 1; 1; 1; 1; 1; 1; 1];
omega_M = [1; 1; 1; 1; 1; 2; 2; 4; 4];
omega_S = [1; 1; 1; 1; 2; 4; 8; 16; 16];

omegas = {omega_C, omega_M, omega_S};

figure;
tiledlayout(3,3);  

for i = 1:3
    omega = omegas{i};

    [t,X] = ode45(@(t,x) SIRVaccineFunction(t, x, beta, gamma, M, N, vax_per_day, h,c,m, omega, true, C_max, dist),[0, T],X0);

    S = X(:,1:n); 
    I = X(:, n+1:2*n);
    R = X(:, 2*n+1:3*n);
    H = X(:, 3*n+1:4*n);
    C = X(:, 4*n+1:5*n);
    D = X(:, 5*n+1:6*n);

    % --- S(t) ---
    nexttile; hold on;
    for k=1:n, plot(t,S(:,k),'LineWidth',1.7,'Color',colors(k,:)); end
    title(['Susceptible - \omega_' num2str(i)]); grid on;

    % --- I(t) ---
    nexttile; hold on;
    for k=1:n, plot(t,I(:,k),'LineWidth',1.7,'Color',colors(k,:)); end
    title(['Infected - \omega_' num2str(i)]); grid on;

    % --- ICU C(t) ---
    nexttile; hold on;
    for k=1:n, plot(t,C(:,k),'LineWidth',1.7,'Color',colors(k,:)); end
    title(['ICU - \omega_' num2str(i)]); grid on;
end

legend(age_groups,'Location','SouthOutside','NumColumns',5);

for k=1:9
    P(k) - (S(end,k) + I(end,k) + R(end,k))
end

%% Con uscite da ospedale e icu
clc

h = [0.001; 0.003; 0.012; 0.032; 0.049; 0.102; 0.166; 0.243; 0.273];
c = [0.050; 0.050; 0.050; 0.050; 0.063; 0.122; 0.274; 0.432; 0.709];
m = [0.00002; 0.00006; 0.00030; 0.00080; 0.00150; 0.00600; 0.02200; 0.05100; 0.09300];

H0 = zeros(n,1);
C0 = zeros(n,1);
M0 = zeros(n,1);
X0 = [S0; I0 ; R0; H0; C0; M0];

[t,X] = ode45(@(t,x)SIRHCMRecoveryFunction(t, x, beta, gamma, M,h,c,m, T_osp, T_icu, N),[0, T],X0);

S = X(:,1:n); 
I = X(:, n+1:2*n);
R = X(:, 2*n+1:3*n);
H = X(:, 3*n+1:4*n);
C = X(:, 4*n+1:5*n);
D = X(:, 5*n+1:6*n);

% --- Grafici H, C, M per fasce ---
figure;
tiledlayout(3,1);

nexttile; hold on;
for k=1:n, plot(t,H(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Hospitalized H(t)'); ylabel('H'); grid on;

nexttile; hold on;
for k=1:n, plot(t,C(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Critical Care C(t)'); ylabel('C'); grid on;

nexttile; hold on;
for k=1:n, plot(t,D(:,k)/P(k),'LineWidth',1.7,'Color',colors(k,:)); end
title('Deaths M(t)'); ylabel('M'); xlabel('Days'); grid on;

legend(age_groups,'Location','SouthOutside','NumColumns',5);

for k=1:9
    P(k) - (S(end,k) + I(end,k) + R(end,k)); % ok
end


%% ======================================================
%       Modello con restrizioni su M in base alle ICU e uscite da ospedale
%       e icu

clc

% NOTA: METTI DEI VALORI SENSATI QUI (E FAI VARIARE dist)
dist = 0.1; % parametro che regola il livello di restrizioni applicate (più è alto, maggiori le restrizioni)
% con 0.1 sembra che ci sia una seconda ondata (HA SENSO CHE SIA MINORE DI
% 1??)
C_max = 14*N/100000; %34.7 * N /100000;

H0 = zeros(n,1);
C0 = zeros(n,1);
M0 = zeros(n,1);
X0 = [S0; I0 ; R0; H0; C0; M0];

[t,X] = ode45(@(t,x)SIRHCMRecoveryFunctionRestriction(t, x, beta, gamma, M,h,c,m, N, dist, C_max, T_osp, T_icu),[0, T],X0);

S = X(:,1:n); 
I = X(:, n+1:2*n);
R = X(:, 2*n+1:3*n);
H = X(:, 3*n+1:4*n);
C = X(:, 4*n+1:5*n);
D = X(:, 5*n+1:6*n);

% --- Grafici restrizioni ---
figure;
tiledlayout(4,1);

nexttile; hold on;
for k=1:n, plot(t,S(:,k)/P(k),'LineWidth',2,'Color',colors(k,:)); end
title('Susceptible with Restrictions'); grid on;

nexttile; hold on;
for k=1:n, plot(t,I(:,k)/P(k),'LineWidth',2,'Color',colors(k,:)); end
title('Infected with Restrictions'); grid on;

nexttile; hold on;
for k=1:n, plot(t,C(:,k)/P(k),'LineWidth',2,'Color',colors(k,:)); end
title('ICU with Restrictions'); xlabel('Days'); grid on;

nexttile; hold on;
plot(t,sum(C,2)/N,'LineWidth',2);
yline(C_max/(dist*N), '--k', 'LineWidth', 2, 'Label', 'C_{max}', 'LabelHorizontalAlignment', 'left');
title('ICU with Restrictions'); xlabel('Days'); grid on;

legend(age_groups,'Location','SouthOutside','NumColumns',5);

for k=1:9
    P(k) - (S(end,k) + I(end,k) + R(end,k)); % ok
end
