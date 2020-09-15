corder = get(gca, 'colororder');
figure(2)
subplot(1, 2, 1)
plot(T, Y2(:,1),'-', 'Color', corder(1, :));
hold on

plot(T, X2(:,1), '-', 'Color', corder(2, :));

plot(T, Yhat2(:,1), '-', 'Color', corder(3, :));

plot(T, YhatIBR2(:, 1), '-', 'Color', corder(4, :));

hold off
xlabel('time')
ylabel('signal value')
title('signal channel 1')


subplot(1, 2, 2)
plot(T, Y2(:, 2), '-', 'Color', corder(1, :));
hold on
plot(T, X2(:, 2), '-', 'Color', corder(2, :));
plot(T, Yhat2(:, 2), '-', 'Color', corder(3, :));
plot(T, YhatIBR2(:, 2), '-', 'Color', corder(4, :));
hold off
xlabel('time')
ylabel('signal value')
title('signal channel 2')
legend('\rho = -0.7 signal', 'observation', '\rho = 0.8 filter', 'IBR filter')