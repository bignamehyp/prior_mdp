function llr = stateToLLR(mu, sigma, policy, prior)
    
     for i = 1 : size(mu,1)
        prob1(i,:) = normcdf(policy(:,2)',mu(i,:),sigma) - normcdf(0, mu(i,:), sigma);
        prob2(i,:) = normcdf(0, mu(i,:),sigma) - normcdf(policy(:,1)',  mu(i,:), sigma);
        llr(i,:) = log(prior / (1 - prior)) + log(prob1(i,:) ./ prob2(i,:));
    end
end