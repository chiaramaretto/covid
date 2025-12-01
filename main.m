clear all
close all
clc

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Simulazione del modello SIR%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model parameters
R_0 = 5.7;
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
% Dati aggregati della popolazione per le 9 fasce (0-9, 10-19, ..., 70-79, 80+)
P = [910147; 963765; 1049452; 1141263; 959404; 951386; 904071; 561572; 265250];
% Calcolo della Popolazione Totale (N)
N = sum(P);

% Calcolo del vettore delle proporzioni p (vettore colonna 9 x 1)
p = P / N;
lambda = max(abs(eigs(M*diag(p))));
gamma = 1/14;           % rate of recovery 
beta = R_0*gamma/lambda;        % rate of infection
T = 350; 
S0 = P;
I0 = ones(n,1);
R0 = zeros(n,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%INSERIRE RISOLUZIONE MODELLO SIR %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0 = [S0; I0 ; R0];

[t,X] = ode45(@(t,x)SIRAgeFunction(t, x, beta, gamma, M, N),[0, T],X0); %con runge-kutta
S = X(:,1:n);
I = X(:, n+1:2*n);
R = X(:, 2*n+1:3*n);
    

% Plots compared with the real data
figure(1)
subplot(2,3,6)
plot(t,I,'r','LineWidth',2);
xlabel('Days'); ylabel('Number of individuals'); title('Infected SIR model')
set(gca,'FontSize',12)

figure(2)
%plot(tt,S,'b',tt,I,'r',tt,R,'g','LineWidth',2);    %eulero
plot(t,S,'b',t,I,'r',t,R,'g','LineWidth',2); 
xlabel('Days'); ylabel('Number of individuals');
legend('S','I','R'); title('SIR model')
set(gca,'FontSize',14)

%%

h = [
    0.001;     % 0-9 anni (0.1%)
    0.003;     % 10-19 anni (0.3%)
    0.012;     % 20-29 anni (1.2%)
    0.032;     % 30-39 anni (3.2%)
    0.049;     % 40-49 anni (4.9%)
    0.102;     % 50-59 anni (10.2%)
    0.166;     % 60-69 anni (16.6%)
    0.243;     % 70-79 anni (24.3%)
    0.273      % 80+ anni (27.3%)
];


c = [
    0.050;     % 0-9 anni (5.0%)
    0.050;     % 10-19 anni (5.0%)
    0.050;     % 20-29 anni (5.0%)
    0.050;     % 30-39 anni (5.0%)
    0.063;     % 40-49 anni (6.3%)
    0.122;     % 50-59 anni (12.2%)
    0.274;     % 60-69 anni (27.4%)
    0.432;     % 70-79 anni (43.2%)
    0.709      % 80+ anni (70.9%)
];


m = [
    0.00002;   % 0-9 anni (0.002%)
    0.00006;   % 10-19 anni (0.006%)
    0.00030;   % 20-29 anni (0.03%)
    0.00080;   % 30-39 anni (0.08%)
    0.00150;   % 40-49 anni (0.15%)
    0.00600;   % 50-59 anni (0.60%)
    0.02200;   % 60-69 anni (2.2%)
    0.05100;   % 70-79 anni (5.1%)
    0.09300    % 80+ anni (9.3%)
];

H0 = zeros(n,1);
C0 = zeros(n,1);
M0 = zeros(n,1);
X0 = [S0; I0 ; R0; H0; C0; M0];


[t,X] = ode45(@(t,x)SIRHCMFunction(t, x, beta, gamma, M,h,c,m, N),[0, T],X0); %con runge-kutta
S = X(:,1:n);
I = X(:, n+1:2*n);
R = X(:, 2*n+1:3*n);
H = X(:, 3*n+1:4*n);
C = X(:, 4*n+1:5*n);
D = X(:, 5*n+1:6*n);

% Plots compared with the real data
figure(1)
plot(t,I,'r','LineWidth',2);
xlabel('Days'); ylabel('Number of individuals'); title('Infected SIR model')
set(gca,'FontSize',12)

figure(2)
%plot(tt,S,'b',tt,I,'r',tt,R,'g','LineWidth',2);    %eulero
plot(t,S,'b',t,I,'r',t,R,'g','LineWidth',2); 
xlabel('Days'); ylabel('Number of individuals');
legend('S','I','R'); title('SIR model')
set(gca,'FontSize',14)

figure(3)
plot(t,H,'b',t,C,'r',t,D,'y',t,R,'g','LineWidth',2); 
xlabel('Days'); ylabel('Number of individuals');
legend('H','C','D', 'R'); title('SIR model')
set(gca,'FontSize',14)