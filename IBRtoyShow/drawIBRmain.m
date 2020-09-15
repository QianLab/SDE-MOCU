clear
% rng default


%parameter setting for SDE
% u = 1;
% theta = (1+u*rand);
% 
% A = 0.01*theta*eye(2);
a = [0;0];
rho = 0.4;
% B = 0.1*[1, rho; rho, 1];
X_initial  = [0; 0];

%blurring filter
TH = 10;
H = 1/TH*ones(1, TH);
H2 = 1;

%noise
sigma  = sqrt(0.01);








%%  sample0 with rho = 0.8
theta = 1;
N = 100;      % # of simulated observations
dt = 1;      % time increment = 1 
rho  = 0.8;
%SDE sample

A = 0.01*theta*eye(2);
B = 0.1*[1, rho; rho, 1];

F = @(t,X) A * X;
G = @(t,X) B;
SDE = sde(F, G, 'StartState', X_initial);
[Y0,T] = SDE.simulate(N, 'DeltaTime', dt, 'nTrials', 1);
X0 = filter(H, H2, Y0)+normrnd(0, sigma, size(Y0));



% plot(T, X0)

%% filter without random parameter and rho = 0.8
rho  = 0.8;
%calculating the variance


[l, r] = meshgrid(0:N, 0:N);
ryx = RYX_factor(r, l, theta, TH);
rxx = RXX_factor(r, l, theta, TH);
RYX = kron([1+rho^2, 2*rho; 2*rho, 1+rho^2], ryx);
RXX = kron([1+rho^2, 2*rho; 2*rho, 1+rho^2], rxx)...
    +sigma^2*kron(eye(2), r==l);
% RYX2 = RYX_matrix_calculate(r, l, theta, rho, TH);
% RXX2 = RXX_matrix_calculate(r, l, theta, rho, TH, sigma);



% G = RYX*(RXX^(-1));
Geight = RYX/RXX; %specific filter
Yhatf0 = Geight*X0(:);
Yhat0 = reshape(Yhatf0, [], 2);
SE_pointeight0 = norm(Yhat0-Y0, 'fro').^2
%% tau filter
rho_t = 0;
RYX_T = kron([1+rho_t^2, 2*rho_t; 2*rho_t, 1+rho_t^2], ryx);
RXX_T = kron([1+rho_t^2, 2*rho_t; 2*rho_t, 1+rho_t^2], rxx)...
    +sigma^2*kron(eye(2), r==l);
G_T = RYX_T/RXX_T;
YTF0 = G_T*X0(:);
YT0 = reshape(YTF0, [], 2);
SE_TF0 = norm(YT0-Y0, 'fro').^2
%% IBR
% M = 10;
% RYXs = zeros(2*length(r), 2*length(l), M);
% RXXs = zeros(2*length(r), 2*length(l), M);
% for m = 1:M
%     tic
% %     theta = 10*rand;
%     rho = 2*rand-1;
% %     rho = rand;
%     RYXs(:, :, m) = RYX_matrix_calculate(r, l, theta, rho, TH);
%     RXXs(:, :, m) = RXX_matrix_calculate(r, l, theta, rho, TH, sigma);    
%     toc
% end
RYXIBR = kron([1+1/3, 0; 0, 1+1/3], ryx);
RXXIBR = kron([1+1/3, 0; 0, 1+1/3], rxx)...
    +sigma^2*kron(eye(2), r==l);
GIBR = RYXIBR/RXXIBR;%IBR filter

YhatIBRf0 = GIBR*X0(:);
YhatIBR0 = reshape(YhatIBRf0, [], 2);
SE_IBR = norm(YhatIBR0-Y0, 'fro').^2


% RYY([2, 3, 4], [15, 56, 2], theta)
% a1 = RYX(1)
% a2 = sqrt(RXX(1)*RYY(1))
% RYX_factor(10, 10, theta, T)

%% sample2 with rho2 = -0.7
rho2 = -0.7;
%SDE sample


B2 = 0.1*[1, rho2; rho2, 1];

F = @(t,X) A * X;
G = @(t,X) B2;
SDE2 = sde(F, G, 'StartState', X_initial);
[Y2,T] = SDE2.simulate(N, 'DeltaTime', dt, 'nTrials', 1);

X2 = filter(H, H2, Y2)+normrnd(0, sigma, size(Y2));

Yhatf2 = Geight*X2(:);
Yhat2 = reshape(Yhatf2, [], 2);
SE_pointeight2 = norm(Yhat2-Y2, 'fro').^2;

YTF2 = G_T*X2(:);
YT2 = reshape(YTF2, [], 2);
SE_TF2 = norm(YT2-Y2, 'fro').^2

YhatIBRf2 = GIBR*X2(:);
YhatIBR2 = reshape(YhatIBRf2, [], 2);
SE_IBR2 = norm(YhatIBR2-Y2, 'fro').^2;


%% draw figure
save toyibr.mat
drawIBR1;
drawIBR2;

%% functions
function RYY = RYY_matrix_calculate(r, l, theta, rho)
    RYY = kron([1+rho^2, 2*rho; 2*rho, 1+rho^2], RYY_factor(r, l, theta));
end

function RYX = RYX_matrix_calculate(r, l, theta, rho, TH)
    RYX = kron([1+rho^2, 2*rho; 2*rho, 1+rho^2], RYX_factor(r, l, theta, TH));
end

function RXX = RXX_matrix_calculate(r, l, theta, rho, TH, sigma)
    RXX = kron([1+rho^2, 2*rho; 2*rho, 1+rho^2], RXX_factor(r, l, theta, TH))...
    +sigma^2*kron(eye(2), r==l);
end


function ryy = RYY_factor(t1, t2, theta)
%return a matrix of size(t1) = size(t2)
%the autocorrelation is calculated by RYY = kron([1+rho^2, 2*rho; 2*rho, 1+rho^2], ryy);
    ryy = 1/(2*theta)*(exp(0.01*theta*(t1+t2))-exp(0.01*theta*abs(t1-t2)));
    ryy = ryy.*(t1>0).*(t2>0);
end



function ryx = RYX_factor(t1, t2, theta, T)
    if (size(t1)~=size(t2))
        error('the size must be the same')
    end
    ryx = zeros(size(t1));
    
    for m = 1:size(t1, 1)
        for n = 1:size(t2, 2)
            fun = @(x) RYY_factor(t1(m, n), x, theta)/T;
            ryx(m, n) = integral(fun, t2(m, n)-T, t2(m, n));
        end
    end
end

function ryx = RYX_factor2(t1, t2, theta, T)
    fun = @(x) RYY_factor(t1, x, theta)/T;
    ryx = integral(fun, t2-T, t2);
end

function rxx = RXX_factor(t1, t2, theta, T)
    if (size(t1)~=size(t2))
        error('the size must be the same')
    end
    rxx = zeros(size(t1));
    fun = @(x, y) RYY_factor(x, y, theta)/(T^2);
    for m = 1:size(t1, 1)
        for n = 1:size(t2, 2)
            rxx(m, n) = integral2(fun, t1(m, n)-T, t1(m, n),...
                t2(m, n)-T, t2(m, n));
        end
    end
end

function rxx = RXX_factor2(t1, t2, theta, T)
    fun = @(x, y) RYY_factor(x, y, theta)/(T^2);
    rxx = integral2(fun, t1-T, t1, t2-T, t2);
end



