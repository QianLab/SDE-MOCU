clear
T = 10;      % duration of simulated observations
dt = 0.01;      % time increment 
rho  = 0;    %the dependence between  winer processes
Y_initial  = [10; 0];
ts0 = cputime;
ts = round(ts0, 3)*1000;
rng(ts);

%SDE parameters
sigma1 = 0.1;
sigma2 = 0.1;
sigmae = 0.2;

% SDE hyper parameters
mu = [0.2;0.5;0.25];
Omega = diag([0.01^2, 0.10^2, 0.02^2]);
fileID = fopen([num2str(ts), 'residual.txt'], 'a');
for m = 1:5

[resk10, resk12, resk21] = IbrResidual(mu, Omega, T, dt, rho, Y_initial, sigma1, sigma2, sigmae, fileID);


fprintf(fileID, '%s\t', 'f');
fprintf(fileID, '%12.12g\t', resk10);
fprintf(fileID, '%12.12g\t', resk12);
fprintf(fileID, '%12.12g\t', resk21);
fprintf(fileID, '\n');

end
fclose(fileID);
% D10 = mean(resk10);
% D12 = mean(resk12);
% D21 = mean(resk21);


% 
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
% 
% 
% % Xt = Yt + et;
% 
% % logS = Y0(2:end, 1, :);
% % logv = Y0(2:end, 2, :);
% % 
% % Y = reshape(logS, 100, []);
% % X = reshape(logv, 100, []);
% % 
% % RYY = (Yt-mean(Yt, 2))*(Yt-mean(Yt, 2))'/MCnum;
% % RYX = (Yt-mean(Yt, 2))*(Xt-mean(Xt, 2))'/MCnum;
% % RXX = (Xt-mean(Xt, 2))*(Xt-mean(Xt, 2))'/MCnum;
% 
% 
% RYY = (Yt)*(Yt)'/MCnum;
% RYX = (Yt)*(Xt)'/MCnum;
% RXX = (Xt)*(Xt)'/MCnum;
% 
% Gtrue = RYX/RXX;%specific filter
% Yhat = Gtrue*Xt(:, 12);
% 
% figure()
% hold on
% for mm = 1:1
% plot( Yt(:, 12), 'r')
% plot( Yhat, 'g')
% % plot(T0, Y0(:, 2, mm))
% % plot(T0, Y0(:, 1, mm)+Y0(:, 2, mm))
% end
% 
% SquareError(Gtrue, Xt, Yt)
% 
% 
% %%
% MM = 1;
% [Y0, T0] = SDE.simulate(T/dt, 'DeltaTime', dt, 'nTrials', MM);
% et0 = normrnd(0, sigmae, time_num, 1);
% X0 = Y0(:, 1)+et0;
% Y0 = Y0(:);
% Yhat0 = Gtrue*X0;
% SquareError(Gtrue, X0, Y0)
% figure()
% hold on
% for mm = 1:MM
% plot( Y0, 'r')
% plot( Yhat0, 'g')
% % plot(T0, Y0(:, 2, mm))
% % plot(T0, Y0(:, 1, mm)+Y0(:, 2, mm))
% end
% % % % % mu = zeros(100, 1);
% % % % % sigma = zeros(100, 1);
% % % % % for m = 1:100
% % % % %     cc = Y0(m, 1, :);
% % % % %     tt = cc(:);
% % % % %     [mu(m), sigma(m)] = normfit(tt);
% % % % % end
% % % % % plot(mu)
% % % % % figure()
% % % % % plot(sigma)
% 
% % t = 0:99;
% % sigma2 =sqrt( 1/200/0.01/theta*(exp(theta*0.01*2*t) - 1)*(1+rho^2));
% % hold on
% % plot(sigma2, 'r')






