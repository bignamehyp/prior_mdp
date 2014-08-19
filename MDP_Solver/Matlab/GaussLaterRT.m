function stopT = GaussLaterRT( d, sigma, mu, nTrials)
stopT = zeros(1,nTrials);
iter = 0;
nT = size(d,1);
thrR = d(:,2)';
thrL = d(:,1)';
while iter < nTrials
    iter = iter + 1;
    o = (randn(1) * sigma + mu);
% 
%     tR = interp1(thrR, 1:nT, o);
%     tL = interp1(thrL, 1:nT, o);
%   
%     stopT(iter) = min(tR,tL);
%     if stopT(iter) == NaN
%         iter = iter - 1;
%     end

    tR  = find(o >= thrR);
    tL = find(o <= thrL);
    if(isempty(tR))
        tR(1) = nT;
    end
    if(isempty(tL))
        tL(1) = nT;
    end
    if(tR(1) > 1)
        rB = tR(1) - 1 + (o - thrR(tR(1)-1)) / (thrR(tR(1)) - thrR(tR(1)-1));
    else
        rB = 1;
    end
    if tL(1) > 1
        lB = tL(1) - 1 + (o - thrL(tL(1)-1)) / (thrL(tL(1)) - thrL(tL(1)-1));
    else
        lB = 1;
    end
    if rB ~= lB
        stopT(iter) = min(rB, lB);
    else
        iter = iter  - 1;
    end
end

end

