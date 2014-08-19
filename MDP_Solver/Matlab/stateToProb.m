function prob = stateToProb(mu, sigma, prior)
    
     for i = 1 : size(mu,1)
        cdf(i,:) = 1 - normcdf(0, mu(i,:), sigma);
        prob(i,:) = prior * cdf(i,:) ./ (prior * cdf(i,:) + (1 - prior) * ( 1 - cdf(i,:)));
    end
end

