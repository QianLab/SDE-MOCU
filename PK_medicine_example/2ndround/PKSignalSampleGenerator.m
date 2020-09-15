function [Xt, Yt] = PKSignalSampleGenerator(McNum, k10, k12, k21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae)

A = [0; 0];
B = [-k12-k10, k21; k12, -k21];
C = [1, rho; rho, 1];

F = @(t,Y) A + B*Y;  %dY = F(Y, t)dt+G(Y, t)dw
G = @(t,Y) [sigma1, 0; 0, sigma2];
SDE = sde(F, G, 'StartState', Y_initial, 'Correlation', C);
[Yt0, t] = SDE.simulate(T/dt, 'DeltaTime', dt, 'nTrials', McNum);
Yt = permute(Yt0, [1, 3, 2]);
Yt = [Yt(:, :, 1);Yt(:, :, 2)];

time_num = size(Yt, 1)/2;
et = normrnd(0, sigmae, time_num, McNum);
Xt = Yt(1:time_num, :)+et;
        
end