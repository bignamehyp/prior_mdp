function [pR RT]= modelRT( d, c, options)
if nargin <= 2
    options = 1;
end


d = d + 1;%matlab index starts from 1
num_c = length(c);
nT = length(d);
tR = zeros(num_c,1);
tL = tR;
pR = tR;
isSampling = zeros(nT,nT + 2);
isRight = zeros(nT,nT + 2);
isLeft = zeros(nT,nT + 2);

for t  = 1:nT
    isSampling(t,d(t,1):d(t,2)) = 1;
    isRight(t,d(t,2)+1:end) = 1;
    isLeft(t, 1:d(t,1)-1) = 1;
end

for i = 1 : num_c
    %u = (c(i) + 1)/2;
    u = (40 * c(i) + 20) / (40 + 20 * c(i));
    filterPDF = zeros(nT, nT + 1); %P(u_{t}, u_{1:t-1} \in sampling)
    condP_R = zeros(1,nT); %P(u_{t+1} \in right, u_{1:t} \in sampling)
    condP_L =  condP_R;  %P(u_{t+1} \in left, u_{1:t} \in sampling)
    filterPDF(1,1) = 1 - u;
    filterPDF(1,2) = u;
    for t = 1 : nT - 1
        filterPDF(t,nT+2) = 0;
        filterPDF(t+1, d(t,1):d(t,2) ) = filterPDF(t+1,d(t,1):d(t,2)) +  filterPDF(t,d(t,1):d(t,2) )* (1-u) .* isSampling(t+1,d(t,1):d(t,2) );
        filterPDF(t+1,(d(t,1):d(t,2) )+1) = filterPDF(t+1,(d(t,1):d(t,2) )+1) + filterPDF(t,d(t,1):d(t,2) ) * u .* isSampling(t+1,(d(t,1):d(t,2))+1);
        condP_R(t+1) = sum(filterPDF(t,d(t,1):d(t,2) ) * (1-u) .* isRight(t+1,d(t,1):d(t,2)) + filterPDF(t,d(t,1):d(t,2) )* u .* isRight(t+1,(d(t,1):d(t,2) )+1));
        condP_L(t+1) = sum(filterPDF(t,d(t,1):d(t,2) )* (1-u) .* isLeft(t+1,d(t,1):d(t,2)) + filterPDF(t,d(t,1):d(t,2) )* u .* isLeft(t+1,(d(t,1):d(t,2) )+1));
    end
    sumR = sum(condP_R);
    sumL = sum(condP_L);
    pR(i) = sumR / (sumR + sumL);
    condP_R = condP_R / sumR;
    condP_L = condP_L / sumL;
    tR(i) = condP_R * (1:nT)';
    tL(i) = condP_L * (1:nT)';
    if options == 1
        tR(i) =  tR(i) / (1 + 0.5 * c(i));
        tL(i) =  tL(i) / (1 + 0.5 * c(i));
    end
end
RT = [tR, tL];
%plot(condP_R(1:2:end));
end

