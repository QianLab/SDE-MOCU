clear
theta = 1;
T = 10;      % duration of simulated observations
dt = 0.01;      % time increment 
rho  = 0;    %the dependence between  winer processes
Y_initial  = [10; 0];


time_num = T/dt+1;
Tarray = 0:dt:T;
% rng(123)

% SDE hyper parameters
mu = [0.2;0.5;0.25];
Omega = diag([0.01^2, 0.1^2, 0.02^2]);


%SDE parameters
sigma1 = 0.1;
sigma2 = 0.1;
sigmae = 0.2;


[IBR_G] = IbrFilterWhole(mu, Omega, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);

% MCnum = 2000;
% A = [0; 0];
% B = [-k12-k10, k21; k12, -k21];
% C = [1, rho; rho, 1];
% 
% F = @(t,Y) A + B*Y;  %dY = F(Y, t)dt+G(Y, t)dw
% G = @(t,Y) [sigma1, 0; 0, sigma2];
% SDE = sde(F, G, 'StartState', Y_initial, 'Correlation', C);
% [Yt0, t] = SDE.simulate(T/dt, 'DeltaTime', dt, 'nTrials', MCnum);
% Yt = permute(Yt0, [1, 3, 2]);
% 
% time_num = size(Yt, 1);
% et = normrnd(0, sigmae, time_num, MCnum);
% 
% Xt = Yt(:, :, 1)+et;
% 
% 
% 
% 
% 
% Yt = [Yt(:, :, 1);Yt(:, :, 2)];


% Xt = Yt + et;

% logS = Y0(2:end, 1, :);
% logv = Y0(2:end, 2, :);
% 
% Y = reshape(logS, 100, []);
% X = reshape(logv, 100, []);
% 


% RYY = (Yt)*(Yt)'/MCnum;
% RYX = (Yt)*(Xt)'/MCnum;
% RXX = (Xt)*(Xt)'/MCnum;
%% filter use
% Gtrue = RYX/RXX;%specific filter
MCnum = 5000;
for m = 1:2
    switch m
        case 1
            K10 = mu(1);
            K12 = mu(2);
            K21 = mu(3);
            [Xt1, Yt1] = PKSignalSampleGenerator(1, K10, K12, K21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
        case 2
            K10 = mu(1)+sqrt(Omega(1, 1))*3;
            K12 = mu(2)+sqrt(Omega(2, 2))*3;
            K21 = mu(3)+sqrt(Omega(3, 3))*3;
            [Xt2, Yt2] = PKSignalSampleGenerator(1, K10, K12, K21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
    end
    [Xt, Yt] = PKSignalSampleGenerator(MCnum, K10, K12, K21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
    RYY = (Yt)*(Yt)'/MCnum;
    RYX = (Yt)*(Xt)'/MCnum;
    RXX = (Xt)*(Xt)'/MCnum;
    switch m
        case 1
            Gmu = RYX/RXX;
            
        case 2
            Gfar = RYX/RXX;
    end
end
