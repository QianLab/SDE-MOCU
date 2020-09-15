function [resk10, resk12, resk21] = IbrResidual(mu, Omega, T, dt, rho, Y_initial, sigma1, sigma2, sigmae, fileID)
%the design value based on mocu
    kMcNum = 50;
    
    %initial
    resk10_array = zeros(kMcNum, 1);
    resk12_array = resk10_array;
    resk21_array = resk10_array;
    
    kSampleSet = mvnrnd(mu, Omega, kMcNum);
    tempn = randi(kMcNum);
    k21 = kSampleSet(tempn, 3);
    for kidx = 1:kMcNum
        disp(kidx)
        K10 = kSampleSet(kidx, 1);
        K12 = kSampleSet(kidx, 2);
        
        [resk10_array(kidx), resk12_array(kidx)] = IbrFilterCost2(K10, K12, k21, mu, Omega, T, dt, rho, Y_initial, sigma1, sigma2, sigmae);

        fprintf(fileID, '%d\t', kidx);
        fprintf(fileID, '%12.12g\t', resk10_array(kidx));
        fprintf(fileID, '%12.12g\t', resk12_array(kidx));
        fprintf(fileID, '%12.12g\t', resk21_array(kidx));
        fprintf(fileID, '\n');
    end
    resk10 = mean(resk10_array);
    resk12 = mean(resk12_array);
    resk21 = mean(resk21_array);
end