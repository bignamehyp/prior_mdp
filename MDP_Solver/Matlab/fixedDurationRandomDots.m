function PR = fixedDurationRandomDots( c, mu, sigma, R_P_H, R_P_L, R_N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 

%Compute the decision boundary, when the state (mu, sigma) satisfies
% mu /sigma > -theta / sigma, the policy choose A_R
% where theta = probit( (R_P_H - R_N_L) / (R_P_H + R_P_L - R_N_H - R_N_L)) 

b = sqrt(2) * erfinv( 2 * (R_P_H - R_N) / ( R_P_H + R_P_L - 2 * R_N) - 1 );

PR = normcdf(c * mu / sigma + b);
 
% numTrials = 100000;
% 
% PR = c;
% 
% %the observation
% for i = 1 : length(c)
%     o = randn(numTrials,1) * sigma + mu * c(i);
%     prob_R = 1.0 - normcdf(0, o, sigma);
%     L_R = prob_R * R_P_H + (1 - prob_R) * R_N;
%     L_L = (1 - prob_R) * R_P_L + prob_R * R_N;
%     PR(i) = 1.0 * sum(L_R >= L_L) / numTrials; 
% end
% 

end

