function [ R_P_min, mu_min, sigma_min, energy_min ] = fitPRRT(PR_EXP, RT_EXP, c, R_P_Set, mu_Set, sigma_Set)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

energy_min = 1.0E18;
weight_PR = 1;


for j = 1 : length(R_P_Set)
    R_P = R_P_Set(j);
    try
        load(sprintf('GaussPolicy_-0.1_%d.0_0.0.txt',R_P));
    catch ME1
        system(sprintf('cd ..; java -jar MDP_Solver.jar -0.1 %d 0; cd Matlab/', R_P));
    end
    d = load(sprintf('GaussPolicy_-0.1_%d.0_0.0.txt',R_P));
    
    for i = 1 : length(mu_Set)
        mu = mu_Set(i);
        for k = 1 : length(sigma_Set)
            sigma = sigma_Set(k);
            
            [PR, RT] = GaussSimulateRT(d,c,mu, sigma);
            
            RT(c < 0,:) = fliplr(RT(c < 0,:));
            linSqrFit_U =  [RT(:,1), ones(length(c),1)] \ RT_EXP;
            RT = RT * linSqrFit_U(1) + linSqrFit_U(2);
            
            R2_PR = double(1 - norm(PR(:,1)./(PR(:,1) + PR(:,2)) - PR_EXP)^2 ...
                / norm(PR_EXP - mean(PR_EXP))^2);
            
            R2_RT = double(1 - norm(RT(:,1) - RT_EXP)^2 ...
                / norm(RT_EXP - mean(RT_EXP))^2);
            
            energy = weight_PR * ( 1 - R2_PR) + (1 - weight_PR) * (1 - R2_RT); 
            sprintf('R_P = %d, mu = %.2f, sigma = %.1f energy = %g\n', R_P, mu,sigma, energy)
            
            if energy < energy_min
                mu_min = mu;
                R_P_min = R_P;
                sigma_min = sigma;
                energy_min = energy;
            end
        end
    end
end



end

