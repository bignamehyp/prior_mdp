clear;
load('monk_data.mat');
c = task(1).coh_set;

d0 = load('policy_-0.1_100_0.txt');
nT = 100;


d = d0(1:nT,1:2);
[PC, RT] = modelRT(d,c);

%[~,~,residue] = regress(task(2).rtc', [RT(:,1) ones(length(c),1)]);
error = log(norm(PC - task(2).pT1', 2) );

tau = 0.1; %temperature
nStart = 100;
nIters = 20;
errorHistory = zeros(nIters,1);
errorHistory(1) = error;

minError = inf;
minD = d;
seqt = [1:nT]';

figure; %Figure for plotting enery history

for iStart = 1: nStart
    randGen = randi(4, nIters * nT,1);
    d = d0(1:nT,1:2);
    error = errorHistory(1);
    for iter = 2 : nIters        
        errorHistory(iter) = 0;
        for i = 1 : nT %Gibbs Sampling
            t = randi(nT-1)+1;
            newD = d;
            currentRand = randGen(t + (iter-1)*nT);
            %Proposed a new move
            if currentRand == 1
                newD(1:t, 1) =  newD(1:t, 1) - 1;
            elseif currentRand == 2
                newD(t:end, 1) =  newD(t:end, 1) + 1;
                
            elseif currentRand == 3
                newD(t:end, 2) =  newD(t:end, 2) - 1;
            else
                newD(1:t, 2) =  newD(1:t, 2) + 1;
            end % end if
            
            newD(newD(:,1) < 0, 1) = 0;
            newD(newD(:,2) > [1:nT]', 1) = seqt(newD(:,2) > [1:nT]');
            
            if sum( newD(:,1) > newD(:,2) ) > 0
                errorHistory(iter) = errorHistory(iter) + error;
                continue;
            end
            
            
            [PC, RT] = modelRT(newD,c);
            %[~,~,residue] = regress(task(2).rtc', [RT(:,1) ones(length(c),1)]);
            %if sum(isnan(residue)) > 0
            %    errorHistory(iter) = errorHistory(iter) + error;
            %    continue;
            %end
            
            newerror = log(norm(PC - task(2).pT1', 2));
            
            
            if newerror < error ||  rand(1) < exp(-(newerror - error)/tau )
                d = newD;
                error = newerror;
            end
            errorHistory(iter) = errorHistory(iter) + error;
            if error < minError
                minError = error;
                minD = d;
            end
            
        end % end of nT
        
        errorHistory(iter) = errorHistory(iter) / nT;
        
        plot(errorHistory(1:iter));
        legend(sprintf('Trial %d',iStart));
        xlim([0,iter]);
        drawnow
    end %MCMC interation
end % different MCMC runs

%Plot the Policy
plotPolicy(minD);
saveas(gcf,'best_fit_policy.fig','fig');
saveas(gcf,'best_fit_policy.jpg','jpg');

%Plot PCRT
[PC, RT] = modelRT(minD,c);
linSqrFit =  [RT(:,1), ones(length(c),1)] \ task(1).rtc';
RT = RT * linSqrFit(1) + linSqrFit(2);

figure;
subplot(2,1,1);
plot(c, PC, '-');
hold on;
plot(c, task(2).pT1,'og',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10);
hold off;
xlim([min(c) max(c)]);
ylabel('Proportional Correct','FontWeight', 'bold','FontSize',30);
set(gca, 'XTick', [-.5 -.25 0 .25 .5],'XMinorTick', 'on',  'LineWidth',2,...
    'FontWeight','bold')
subplot(2,1,2);
%Show reaction times of correct trials only
plot(c, RT(:,1),'-');
hold on;
errorbar(c, task(2).rtc, task(2).rtc_se, 'go', ...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10);hold off;
xlim([min(c) max(c)]);
xlabel('Motion Strength','FontWeight', 'bold','FontSize',30);
ylabel('Reaction Time','FontWeight', 'bold','FontSize',30);
set(gca, 'XTick', [-.5 -.25 0 .25 .5],'XMinorTick', 'on', 'LineWidth',2,...
    'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[12 15]);
set(gcf,'paperposition',[0,0,12,15]);
saveas(gcf,'best_fit_prior.fig','fig');
saveas(gcf,'best_fit_prior.jpg','jpg');
save bestfit.mat
