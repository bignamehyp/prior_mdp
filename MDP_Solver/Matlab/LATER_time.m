function stopT = LATER_time( d, tau )
nTrials = 10000;
stopT = zeros(1,nTrials);
iter = 0;
nT = size(d,1);
thr = d(:,2) ./ (1:nT)';
%prob = 0.5 - 0.5 * exp(-(1:nT)/tau);
prob = 0.5*ones(1,nT);
prob(1:tau) = 0;
while iter < nTrials
    iter = iter + 1;
    o = rand(1,nT);
    o(o > prob)  = 1;
    o(o <= prob) = 0;
    
    o = cumsum(o) ./ (1:nT);
    temp  = find(o > thr');
    if(~isempty(temp))
        stopT(iter) = temp(1);
    else
        iter = iter  - 1;
    end
end
    bins = (nT-10):-8:4;
    hist(1.0./stopT,1.0./bins);
end

