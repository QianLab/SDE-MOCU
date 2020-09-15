function [resk10, resk12, resk21] = IbrResidual(mu, Omega, T, dt, rho, Y_initial, sigma1, sigma2, sigmae)
%the design value based on mocu
    kMcNum = 200;
    y1only = 0;
    cost_direct = 1;
    
    %initial
    resk10_array = zeros(kMcNum, 1);
    resk12_array = resk10_array;
    resk21_array = resk10_array;
    
    kSampleSet = mvnrnd(mu, Omega, kMcNum);
    fileID = fopen('residual-record.txt', 'a');
    for kidx = 1:kMcNum
        disp(kidx)
        K10 = kSampleSet(kidx, 1);
        K12 = kSampleSet(kidx, 2);
        K21 = kSampleSet(kidx, 3);
        
        if cost_direct == 1
            [resk10_array(kidx), resk12_array(kidx), resk21_array(kidx)] = IbrFilterCost(K10, K12, K21, mu, Omega, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
        else
        
            [Xt, Yt] = PKSignalSampleGenerator(5000, K10, K12, K21, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);
            [IBR_G_k10, IBR_G_k12, IBR_G_k21] = IbrFilter(K10, K12, K21, mu, Omega, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);

            if y1only == 0
                resk10_array(kidx) = SquareError(IBR_G_k10, Xt, Yt);%if we have to calculate IBR filter with large computation, it should be better to sample more x and y
                resk12_array(kidx) = SquareError(IBR_G_k12, Xt, Yt);
                resk21_array(kidx) = SquareError(IBR_G_k21, Xt, Yt);
            else
                resk10_array(kidx) = SquareError(IBR_G_k10(1:1001, :), Xt, Yt(1:1001, :));%if we have to calculate IBR filter with large computation, it should be better to sample more x and y
                resk12_array(kidx) = SquareError(IBR_G_k12(1:1001, :), Xt, Yt(1:1001, :));
                resk21_array(kidx) = SquareError(IBR_G_k21(1:1001, :), Xt, Yt(1:1001, :));
            end
        end
        
        fprintf(fileID, '%d\t', kidx);
        fprintf(fileID, '%12.12g\t', resk10_array(kidx));
        fprintf(fileID, '%12.12g\t', resk12_array(kidx));
        fprintf(fileID, '%12.12g\t', resk21_array(kidx));
        fprintf(fileID, '\n');
    end
    fclose(fileID);
    resk10 = mean(resk10_array);
    resk12 = mean(resk12_array);
    resk21 = mean(resk21_array);
end