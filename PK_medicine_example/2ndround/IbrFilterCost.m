function [k10cost, k12cost, k21cost] = IbrFilterCost(k10, k12, k21, mu, Omega, T, dt, rho, Y_initial, sigma1, sigma2, sigmae)
%a PK model ibr filter, calculated in a relatively general Monte Carlo way

%we should sample posterior given k10 in general, but here k10 k12 k21 are
%independent
kMcNum = 3000;
SNum = 50;
kSampleSet = mvnrnd(mu, Omega, kMcNum);

Nx = T/dt+1;
Ny = 2*Nx;

%initial
% rxx_G_k10_sum = zeros(Nx, Nx);
% ryx_G_k10_sum = zeros(Ny, Nx);
Xtcon = zeros(Nx, kMcNum*SNum);
Ytcon = zeros(Ny, kMcNum*SNum);

for m = 1:3
    for kidx = 1:kMcNum
        K10 = kSampleSet(kidx, 1);
        K12 = kSampleSet(kidx, 2);
        K21 = kSampleSet(kidx, 3);
        switch m
            case 1
                [Xt, Yt] = PKSignalSampleGenerator(SNum, k10, K12, K21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
            case 2
                [Xt, Yt] = PKSignalSampleGenerator(SNum, K10, k12, K21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
            case 3
                [Xt, Yt] = PKSignalSampleGenerator(SNum, K10, K12, k21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
        end
        Xtcon(:, kidx*SNum-SNum+1:kidx*SNum) = Xt;
        Ytcon(:, kidx*SNum-SNum+1:kidx*SNum) = Yt;
    end
    switch m
        case 1
            [k10cost, IBR_k10] = IBRCostcalc(Xtcon, Ytcon);
        case 2
            [k12cost, IBR_k12] = IBRCostcalc(Xtcon, Ytcon);
        case 3
            [k21cost, IBR_k21] = IBRCostcalc(Xtcon, Ytcon);
    end    
end

for m = 1:3
    for kidx = 1:kMcNum
        K10 = kSampleSet(kidx, 1);
        K12 = kSampleSet(kidx, 2);
        K21 = kSampleSet(kidx, 3);
        switch m
            case 1
                [Xt, Yt] = PKSignalSampleGenerator(SNum, k10, K12, K21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
            case 2
                [Xt, Yt] = PKSignalSampleGenerator(SNum, K10, k12, K21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
            case 3
                [Xt, Yt] = PKSignalSampleGenerator(SNum, K10, K12, k21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
        end
        Xtcon(:, kidx*SNum-SNum+1:kidx*SNum) = Xt;
        Ytcon(:, kidx*SNum-SNum+1:kidx*SNum) = Yt;
    end
    switch m
        case 1
            SquareError(IBR_k10, Xtcon, Ytcon)
        case 2
            SquareError(IBR_k12, Xtcon, Ytcon)
        case 3
            SquareError(IBR_k21, Xtcon, Ytcon)
    end
end



end

function [ibr_cost, IBR_G_k] = IBRCostcalc(Xt, Yt)
%     ryx_G_k_sum = (Yt)*(Xt)';
%     rxx_G_k_sum = (Xt)*(Xt)';
%     ryy_G_k_sum = (Yt)*(Yt)';
    ryxIBR_G_k = (Yt)*(Xt)'./size(Xt, 2);
    rxxIBR_G_k = (Xt)*(Xt)'./size(Xt, 2);
    ryyIBR_G_k = (Yt)*(Yt)'./size(Xt, 2);
     IBR_G_k = ryxIBR_G_k/rxxIBR_G_k;
    ibr_cost = trace(ryyIBR_G_k) - trace((ryxIBR_G_k/rxxIBR_G_k)*ryxIBR_G_k');
end