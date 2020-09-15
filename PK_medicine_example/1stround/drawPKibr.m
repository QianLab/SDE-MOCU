%% draw figure
load('PKibr')
corder = get(gca, 'colororder');
Gtau = Gmu;
for m = 1:2
    figure(m)
    hold on
    switch m
        case 1
            Xt0 = Xt1;
            Yt0 = Yt1;
            signalstring = '$\theta = \mu$ signal';
        case 2
            Xt0 = Xt2;
            Yt0 = Yt2;
            signalstring = '$\theta = \mu+3\sigma_\theta$ signal';
    end
    Yhat1 = Gfar*Xt0;
    Yhat2 = Gtau*Xt0;
%     Yhat2 = Gfar*Xt0;
    Yhat3 = IBR_G*Xt0;
    SquareError(Gmu, Xt0, Yt0)
    SquareError(Gfar, Xt0, Yt0)
    SquareError(IBR_G, Xt0, Yt0)
    for n = 1:2
        switch n
            case 1
                subplot(1, 2, 1)
                title('central compartment', 'FontSize', 15)
                idx = 1:time_num;
            case 2
                subplot(1, 2, 2)
                title('peripheral compartment', 'FontSize', 15)
                idx = time_num+1:2*time_num;
        end
        hold on
        plot(Tarray,  Yt0(idx), 'Color', corder(1, :), 'LineWidth', 2)
        plot(Tarray, Yhat1(idx),  'Color', corder(2, :), 'LineWidth', 2) 
        plot(Tarray, Yhat2(idx), 'Color', corder(5, :), 'LineWidth', 2)   
        plot(Tarray, Yhat3(idx), 'Color', corder(3, :), 'LineWidth', 2)   
        xlabel('time', 'FontSize', 15)
        ylabel('signal value', 'FontSize', 15)
        
        
        set(gcf, 'Position',  [100, 100, 800, 380])
        if n == 1
            plot(Tarray(1:20:end), Xt0(1:20:end), '+', 'Color', corder(4, :))%, 'MarkerSize',2) 
            leg = legend(signalstring, '$\theta = \mu+3\sigma_\theta$ filter', '$\tau$-robust filter', 'IBR filter','observation', 'FontSize', 15)
            set(leg, 'Interpreter', 'latex')
        end
        hold off
    end
end
    
% Yhat = IBR_G*Xt;
% 
% figure()
% hold on
% for mm = 1:1
% plot( Yt, 'r')
% plot( Yhat, 'g')
% % plot(T0, Y0(:, 2, mm))
% % plot(T0, Y0(:, 1, mm)+Y0(:, 2, mm))
% end
% 
% % SquareError(Gtrue, Xt, Yt)


%%
% % % 
% % % kSampleSet = mvnrnd(mu, Omega, 1);
% % % % k10 = kSampleSet(1);
% % % % k12 = kSampleSet(2);
% % % % k21 = kSampleSet(3);
% % % k10 = 0.22;
% % % k12 = 0.5;
% % % k21 = 0.29;
% % % 
% % % 
% % % MCnum = 20;
% % % A = [0; 0];
% % % B = [-k12-k10, k21; k12, -k21];
% % % C = [1, rho; rho, 1];
% % % 
% % % F = @(t,Y) A + B*Y;  %dY = F(Y, t)dt+G(Y, t)dw
% % % G = @(t,Y) [sigma1, 0; 0, sigma2];
% % % SDE = sde(F, G, 'StartState', Y_initial, 'Correlation', C);
% % % [Yt0, t] = SDE.simulate(T/dt, 'DeltaTime', dt, 'nTrials', MCnum);
% % % Yt = permute(Yt0, [1, 3, 2]);
% % % 
% % % time_num = size(Yt, 1);
% % % et = normrnd(0, sigmae, time_num, MCnum);
% % % 
% % % Xt = Yt(:, :, 1)+et;
% % % 
% % % 
% % % 
% % % 
% % % 
% % % Yt = [Yt(:, :, 1);Yt(:, :, 2)];
% % % 
% % % 
% % % % % Xt = Yt + et;
% % % % 
% % % % logS = Y0(2:end, 1, :);
% % % % logv = Y0(2:end, 2, :);
% % % % 
% % % % Y = reshape(logS, 100, []);
% % % % X = reshape(logv, 100, []);
% % % 
% % % 
% % % [Xt, Yt] = PKSignalSampleGenerator(MCnum, k10, k12, k21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
% % % figure()
% % % plot(Yt)
% % % 
% % % 
% % % RYY = (Yt)*(Yt)'/MCnum;
% % % RYX = (Yt)*(Xt)'/MCnum;
% % % RXX = (Xt)*(Xt)'/MCnum;
% % % Gtrue = RYX/RXX;%specific filter
% % % MM = 1;
% % % [Y0, T0] = SDE.simulate(T/dt, 'DeltaTime', dt, 'nTrials', MM);
% % % et0 = normrnd(0, sigmae, time_num, 1);
% % % X0 = Y0(:, 1)+et0;
% % % Y0 = Y0(:);
% % % Yhat0 = Gtrue*X0;
% % % SquareError(Gtrue, X0, Y0)
% % % figure()
% % % hold on
% % % for mm = 1:MM
% % % plot( Y0, 'r')
% % % plot( Yhat0, 'g')
% % % end