function [ muMin prMin, sigmaMin, energy ] = fitBiasData(task )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%[mu sigma, PR] = (1,  0.9, 80) for ssp
%               = (1.5, 1.1, 100) for monk
%               = (1.2 0.9, 80) for Lsp
%[mu, P_R]      =  (1.5  180) for Sacc
%[mu, P_R]      =  (1.5, 220) for Lacc


muSet = 0.9:0.1:1.2;
PrSet = 60:20:100;
sigmaSet = 0.6:0.1:1.5;
c = task(1).coh_set;
energy = zeros(length(PrSet), length(muSet));
minE = 1.0^8;

for j = 1 : length(PrSet)
    P_R = PrSet(j);
    try
        load(sprintf('GaussPolicy_-0.1_%d.0_0.0.txt',P_R));
    catch ME1
        system(sprintf('cd ..; java -jar MDP_Solver.jar -0.1 %d 0; cd Matlab/', P_R));
    end
    d0_U = load(sprintf('GaussPolicy_-0.1_%d.0_0.0.txt',P_R));
    for i = 1 : length(muSet)
        mu = muSet(i);
        for k = 1 : length(sigmaSet)
            sigma = sigmaSet(k);
            [PR_U, Steps_U] = GaussSimulateRT(d0_U,c,mu, sigma);
            Steps_U(c < 0,:) = fliplr(Steps_U(c < 0,:));
            linSqrFit_U =  [Steps_U(:,1), ones(length(c),1)] \ task(1).rtc';
            RT_U = Steps_U * linSqrFit_U(1) + linSqrFit_U(2);
            
            R2PRU = double(1 - norm(PR_U(:,1)./(PR_U(:,1) + PR_U(:,2)) - task(1).pT1')^2 / norm(task(1).pT1' - mean(task(1).pT1))^2);
            R2TimeU = double(1 - norm( (RT_U(:,1) - task(1).rtc'))^2 / norm(task(1).rtc' - mean(task(1).rtc))^2);
            
            %energy(i,j) =  1 - R2PRU + 1 - R2TimeU;
            energy(i,j) = 1 - R2PRU;
            sprintf('R_P = %d, mu = %.2f, sigma = %.1f energy = %g\n', P_R, mu,sigma, energy(i,j))
            if energy(i,j) < minE
                muMin = mu;
                prMin = P_R;
                sigmaMin = sigma;
                minE = energy(i,j);
            end
        end
    end
end
end

