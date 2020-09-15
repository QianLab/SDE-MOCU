
% Theta             Rho         mocurho        mocutheta             
A = [    0.500000,     0.500000,     1.198307,     1.119499
    0.500000,     1.500000,     1.060496,     1.038684
    0.500000,     5.000000,     0.966690,     0.983833
    1.500000,     0.500000,     1.219811,     1.153601
    1.500000,     1.500000,     1.061642,     1.058926
    1.500000,     5.000000,     0.961472,     1.000380
    2.000000,     0.500000,     1.201162,     1.150782
    2.000000,     1.500000,     1.069263,     1.071604
    2.000000,     5.000000,     0.967818,     1.012086];
thetapara = reshape(A(:, 1), 3, 3);
thetavar = thetapara.^2/12;
rhopara = reshape(A(:, 2), 3, 3);
rhovar = 1./(2*rhopara);
mocurho = reshape(A(:, 3),3 , 3);
mocutheta = reshape(A(:, 4), 3, 3);

[xq, yq] = meshgrid(linspace(thetavar(1), thetavar(end), 20), linspace(rhovar(1), rhovar(end), 20));
rhoq = griddata(thetavar, rhovar, mocurho, xq, yq, 'cubic');
figure()
thetaq = griddata(thetavar, rhovar, mocutheta, xq, yq, 'cubic');
colormap([cool(64); hot(64)])
s1 = surf(xq, yq, rhoq, 'FaceColor', 'interp', 'EdgeColor', 'none');

hold on
s2 = surf(xq, yq, thetaq, 'FaceColor', 'interp','EdgeColor', 'none')
plot3(thetavar, rhovar, mocurho, 'bo')
plot3(thetavar, rhovar, mocutheta, 'ro')
hold off
m = 64;
cmin = min(rhoq(:));
cmax = max(rhoq(:));
C1 = min(m, round((m-1)*(rhoq-cmin)/(cmax-cmin))+1);
C2 = 64+C1;

set(s1, 'CData', C1);
set(s2, 'CData', C2);
caxis([min(C1(:)) max(C2(:))])
colorbar('eastoutside','Ticks', [32, 96], 'TickLabels', {'R(\theta_{1})', 'R(\theta_2)'}, 'FontSize', 15)
% c.Label.String = 'D(\rho)                 D(\theta)'
xlabel('Var(\theta_1)', 'FontSize', 15)
ylabel('Var(\theta_2)', 'FontSize', 15)
zlabel('R(\theta_1) or R(\theta_2)', 'FontSize', 15)


