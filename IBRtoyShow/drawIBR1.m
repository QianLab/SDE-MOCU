load('toyibr.mat')
corder = get(gca, 'colororder');
for n = 1:2
figure(n)
switch n
    case 1
        Y = Y0;
        Yhat = Yhat0;
        YhatIBR = YhatIBR0;
        X = X0;
        YT = YT0;
        
    case 2
        Y = Y2;
        Yhat = Yhat2;
        YhatIBR = YhatIBR2;
        YT = YT2;
        X = X2;
end
for m = 1:2
subplot(1, 2, m)
plot(T, Y(:,m),'-', 'Color', corder(1, :), 'LineWidth', 2);
hold on



plot(T, Yhat(:,m), '-', 'Color', corder(2, :), 'LineWidth', 2);
plot(T, YT(:, m), '-', 'Color', corder(5, :), 'LineWidth', 2);
plot(T, YhatIBR(:, m), '-', 'Color', corder(3, :), 'LineWidth', 2);
plot(T, X(:,m), '+', 'Color', corder(4, :));


hold off
xlabel('time', 'FontSize', 15)
% xlim([40, 60])
switch m
    case 1
        title('signal channel 1','FontSize', 15)
        ylabel('signal value','FontSize', 15)
    case 2
        title('signal channel 2','FontSize', 15)
        switch n
            case 1
                leg = legend('$\theta_2$ = 0.8 signal',  '$\theta_2$ = 0.8 filter','$\tau$-filter', 'IBR filter', 'observation','FontSize', 12);

            case 2
                leg = legend('$\theta_2$ = -0.7 signal',  '$\theta_2$ = 0.8 filter','$\tau$-filter', 'IBR filter', 'observation', 'FontSize', 12);
        end
        set(leg, 'Interpreter', 'latex')
end

end
set(gcf, 'Position',  [100, 100, 800, 380])
end


% subplot(1, 2, 2)
% plot(T, Y0(:, 2), '-', 'Color', corder(1, :));
% hold on
% plot(T, X0(:, 2), '-', 'Color', corder(2, :));
% plot(T, Yhat0(:, 2), '-', 'Color', corder(3, :));
% plot(T, YhatIBR0(:, 2), '-', 'Color', corder(4, :));
% hold off
% xlabel('time')
% ylabel('signal value')
% title('signal channel 2')
% 