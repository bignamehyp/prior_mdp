function bias = dynamicBias(isSave)
if nargin <= 1
    isSave = 0;
end

d0_U = load('GaussPolicy_-0.1_100.0_0.0.txt');

priorSet = [0.6 0.7 0.8 0.9];
nT = size(d0_U,1);
muSet = -2:0.01:2;
m = length(muSet);

bias = zeros(length(priorSet),nT);
legends = [];
for p = 1: length(priorSet)
    prior = priorSet(p);
    disp(prior);
    legends{p}=sprintf('Prior = %.1f', prior);
    filename = sprintf('GaussPolicy_biased_%.1f.txt', prior);
    d0_B = load(filename);
    prB = zeros(m,nT);
    prU = zeros(m,nT);

    prB(muSet >= d0_B(nT,2), nT) = 1;
    prU(muSet >= d0_U(nT,2), nT) = 1;
    
    bias(p,nT) = log(sum(prB(:,nT)) / (  m - sum(prB(:,nT)))) / log( sum(prU(:,nT)) / ( m - sum(prU(:,nT))) );

    for t =  nT - 1: -1 : 1
        
        for i = 1: m
            mu = muSet(i);
            
            if mu >= d0_U(t,2)
                prU(i,t) = 1;
            elseif mu > d0_U(t,1)
                tran = normpdf(muSet, mu, 1.0/t + 1);
                tran = tran / sum(tran);
                prU(i,t) =  sum(tran .* prU(:, t+1)');
            else
                prU(i,t) = 0;
            end
                        
            if mu >= d0_B(t,2)
                prB(i,t) = 1;
            elseif mu > d0_B(t,1)
                tran = normpdf(muSet,mu, 1.0/t+1) ...
                    .* (prior + (1 - prior * 2) * normcdf(0, (mu * t + muSet)/(t+1), 1.0/(t+1))) ...
                    / (prior + ( 1- prior * 2) * normcdf(0, mu, 1.0/t));
                tran = tran/ sum(tran);
                prB(i,t) = sum( tran .* prB(:,t+1)') ;
            else
                prB(i,t)= 0;
            end
        end
        
        xSet  = 0.1:0.1:1;
        sumB = 0;
        sumU = 0;
        sumB2 = 0;
        sumU2 = 0;
        for x = xSet
            sumB =  sumB + sum(prB(:,t) .* normpdf(muSet, x,  1.0/t )')/sum(normpdf(muSet, x, 1.0/t ));
            sumB2 =  sumB2 + sum((1 - prB(:,t)) .* normpdf(muSet, x,  1.0/t )')/sum(normpdf(muSet, x, 1.0/t ));
            sumU =  sumU + sum(prU(:,t) .* normpdf(muSet, x, 1.0/t )')/sum(normpdf(muSet, x,  1.0/t ));
            sumU2 =  sumU2 + sum((1 - prU(:,t)) .* normpdf(muSet, x,  1.0/t )')/sum(normpdf(muSet, x, 1.0/t ));
        end
        sumB = sumB/(sumB + sumB2);
        sumU = sumU/(sumU + sumU2);

        %bias(p,t) =  log(sumB / (1 - sumB)) / log( sumU / (1 - sumU) );
        bias(p,t) = sumU;
    end
end


figure;
plot(bias');
xlim([1 nT]);
legend(legends,'Location','Best');
xlabel('Time Step t','FontWeight', 'bold','FontSize',30);
ylabel('Bias Signal','FontWeight', 'bold','FontSize',30);
set(gca, 'LineWidth',2,'FontSize',30, 'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[8 8]);
set(gcf,'paperposition',[0,0,8,8]);
if isSave
    saveas(gcf,'DynamicBias.fig','fig');
    saveas(gcf,'DynamicBias.jpg','jpg');
end

end