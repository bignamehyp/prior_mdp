%% Fix Duration LIP Activities
% Decision boundaries when LIP 

%% Reaction Time LIP Activities
prior = 0.9;
d0_B = load('Policy/GaussPolicy_0.90_-0.1_100.0_0.0_30.txt');
d0_B2 = load('Policy/GaussPolicy_0.10_-0.1_100.0_0.0_30.txt');

d0_U = load('Policy/GaussPolicy_0.50_-0.1_100.0_0.0_30.txt');
%d0_U2 = load('Policy/GaussPolicy_0.50_-0.1_50.0_0.0_30.txt');

nT = size(d0_B,1);
 
mu = 1.4;
sigma = 1;
c = 0.256;
 

boundary_HH = stateToProb(d0_U(1:nT,2)', sigma./sqrt(1:nT), 0.5);
boundary_HH = log10(boundary_HH ./(1 - boundary_HH));

% boundary_LL = stateToProb(d0_U2(1:nT,2)', sigma./sqrt(1:nT), 0.5);
% boundary_LL = log10(boundary_LL ./(1 - boundary_LL));

boundary_HL = stateToProb(d0_B(1:nT,2)', sigma./sqrt(1:nT), prior);
boundary_HL = log10(boundary_HL ./ ( 1 - boundary_HL));
 
boundary_LH = stateToProb(d0_B2(1:nT,2)', sigma./sqrt(1:nT), 1 - prior);
boundary_LH = log10(boundary_LH ./ ( 1 - boundary_LH));


[~, ~, ~, state_seqs_HH] = GaussSimulateRT( d0_U(1:nT,:), c, mu, sigma);
LIP_HH = stateToProb(state_seqs_HH, sigma./sqrt(1:nT), 0.5);
LIP_HH = log10(LIP_HH ./ ( 1 - LIP_HH));
LIP_HH =   LIP_HH - repmat(boundary_HH, length(c),1);

% [~, ~, ~, state_seqs_LL] = GaussSimulateRT( d0_U2(1:nT,:), c, mu, sigma);
% LIP_LL = stateToProb(state_seqs_LL, sigma./sqrt(1:nT), 0.5);
% LIP_LL = log10(LIP_LL ./ ( 1 - LIP_LL));
% LIP_LL =   LIP_LL - repmat(boundary_LL, length(c),1);

[~, ~, ~, state_seqs_HL] = GaussSimulateRT( d0_B(1:nT,:), c, mu, sigma,prior);

LIP_HL = stateToProb(state_seqs_HL, sigma./sqrt(1:nT), prior);
LIP_HL = log10(LIP_HL ./ (1 - LIP_HL));
LIP_HL =   LIP_HL - repmat(boundary_HL, length(c),1);
 
[~, ~, ~, state_seqs_LH] = GaussSimulateRT( d0_B2(1:nT,:), c, mu, sigma, 1- prior);

LIP_LH = stateToProb(state_seqs_LH, sigma./sqrt(1:nT),  1- prior);
LIP_LH = log10(LIP_LH ./ (1 - LIP_LH));
LIP_LH =   LIP_LH - repmat(boundary_LH, length(c),1);

rangeLIP = -min(LIP_LH);

LIP_HH = (LIP_HH + rangeLIP)/rangeLIP;
LIP_LH = (LIP_LH + rangeLIP)/rangeLIP;
LIP_HL = (LIP_HL + rangeLIP)/rangeLIP;


 



%% Reaction Time LIP Activities
prior = 0.7;
d0_B = load('Policy/GaussPolicy_0.70_-0.1_100.0_0.0_100.txt');
d0_U = load('Policy/GaussPolicy_0.50_-0.1_100.0_0.0_100.txt');
d0_B2 = d0_B;
d0_B2(:,1) = -d0_B(:,2);
d0_B2(:,2) = -d0_B(:,1);

nT = size(d0_B,1);
 
mu = 1;
sigma = 1;
c = [0.032 0.064 0.128 0.256];
 

boundary_U = stateToProb(d0_U(1:nT,2)', sigma./sqrt(1:nT), 0.5);
boundary_U = log10(boundary_U ./(1 - boundary_U));
 
boundary_B = stateToProb(d0_B(1:nT,2)', sigma./sqrt(1:nT), prior);
boundary_B = log10(boundary_B ./ ( 1 - boundary_B));
 
boundary_B2 = stateToProb(d0_B2(1:nT,2)', sigma./sqrt(1:nT), 1-prior);
boundary_B2 = log10(boundary_B2 ./ ( 1 - boundary_B2));

% boundary_B =  boundary_U - 0.4 * (1 - exp(-(1:nT)/nT*16)); 
% cdf = 1.0 ./ (1.0 + 10.^(boundary_B- log10(prior / ( 1 - prior) ))); 
% cdf = erfinv(1 - 2 * cdf) * sqrt(2) * sigma./sqrt(1:nT);
% d0_B(1:nT,1) = d0_B(1:nT,1) + d0_B(1:nT,2) - cdf';
% d0_B(1:nT,2) = cdf;
%  

[~, ~, ~, state_seqs_U] = GaussSimulateRT( d0_U(1:nT,:), c, mu, sigma);
LIP_U = stateToProb(state_seqs_U, sigma./sqrt(1:nT), 0.5);
LIP_U = log10(LIP_U ./ ( 1 - LIP_U));
[~, spanU] = max(LIP_U,[],2);
LIP_U =   LIP_U - repmat(boundary_U, length(c),1);

[~, ~, ~, state_seqs_B] = GaussSimulateRT( d0_B(1:nT,:), c, mu, sigma,prior);

LIP_B = stateToProb(state_seqs_B, sigma./sqrt(1:nT), prior);
LIP_B = log10(LIP_B ./ (1 - LIP_B));
[~, spanB] = max(LIP_B,[],2);
LIP_B =   LIP_B - repmat(boundary_B, length(c),1);

[~, ~, ~, state_seqs_B2] = GaussSimulateRT( d0_B2(1:nT,:), c, mu, sigma,1-prior);

LIP_B2 = stateToProb(state_seqs_B2, sigma./sqrt(1:nT), 1 - prior);
LIP_B2 = log10(LIP_B2 ./ (1 - LIP_B2));
[~, spanB2] = max(LIP_B2,[],2);
LIP_B2 =   LIP_B2 - repmat(boundary_B2, length(c),1);

rangeLIP = -min(min(LIP_U));

LIP_U = (LIP_U + rangeLIP) / rangeLIP;
LIP_B = (LIP_B + rangeLIP) / rangeLIP;
LIP_B2 = (LIP_B2 + rangeLIP) / rangeLIP;

if size(d0_B,1)  >= 100
    for i = 1 : length(spanU)
        LIP_U(i, spanU(i)+1 : nT) = nan;
        LIP_B(i, spanB(i)+1 : nT) = nan;
        LIP_B2(i, spanB2(i)+1 : nT) = nan;

    end     
end



%  
% buildup_B = zeros(1,length(c));
% buildup_U = zeros(1,length(c));
% 
% 
% %Build up Rate
% for i = 1 : length(c) 
%     x = 1:spanB(i);
%     y = LIP_B(i,x)  ;
%     buildup_B(i) =   (mean( x .* y ) - mean(x) * mean(y)) / (mean(x.^2) - mean(x).^2);
%     %x = 1:spanU(i);
%     y = LIP_U(i,x);
%     buildup_U(i) =  (mean( x .* y ) - mean(x) * mean(y)) / (mean(x.^2) - mean(x).^2);
% end

% subplot(1,2,2);
% plot(c, buildup_B,'r','MarkerSize',10,'Marker','o','DisplayName', '67:33 prior');
% hold on;
% plot(c, buildup_U,'b','MarkerSize',10,'Marker','o','DisplayName', '50:50 prior');
% xlabel('Coherence');
% ylabel('Buildup rate');
% xlim([c(1), c(end)]);
% legend('show','Location','NorthWest');
% set(gcf,'paperunits','inches');
% set(gcf,'papersize',[9 15]);
% set(gcf,'paperposition',[0,0,9,15]);
 
 
%%

filename = 'monk';
load(sprintf('%s_data.mat',filename));
c = task(1).coh_set;

mu = 1;
sigma = 1;

boundary_B2 =  boundary_U - 0.12 +  0.15*(1 :nT).^1.4/nT^1.4; 
cdf = 1.0 ./ (1.0 + 10.^(boundary_B2 - log10(prior / ( 1 - prior) ))); 
cdf = erfinv(1 - 2 * cdf) * sqrt(2) * sigma./sqrt(1:nT);
d0_B(1:nT,1) = d0_B(1:nT,1) + d0_B(1:nT,2) - cdf';
d0_B(1:nT,2) = cdf;

[PR_U, RT_HL] = GaussSimulateRT(d0_U,c,mu, sigma);
[PR_B, RT_LH] = GaussSimulateRT(d0_B,c,mu, sigma);

RT_HL(c < 0,:) = fliplr(RT_HL(c < 0,:));
RT_LH(c < 0,:) = fliplr(RT_LH(c < 0,:));

linSqrFit_U =  [RT_HL(:,1), ones(length(c),1)] \ task(1).rtc';

linSqrFit_B =  [RT_LH(:,1), ones(length(c),1)] \ task(2).rtc';

RT_U = RT_HL * linSqrFit_U(1) + linSqrFit_U(2);
RT_B = RT_LH * linSqrFit_U(1) + mean(task(2).rtc' - RT_LH(:,1) * linSqrFit_U(1));
%RT_B = RT_LH * linSqrFit_B(1) + linSqrFit_B(2);

R2PRU = double(1 - norm(PR_U(:,1)./(PR_U(:,1) + PR_U(:,2)) - task(1).pT1')^2 / norm(task(1).pT1' - mean(task(1).pT1))^2);
R2TimeU = double(1 - norm( (RT_U(:,1) - task(1).rtc'))^2 / norm(task(1).rtc' - mean(task(1).rtc))^2);
R2PRB = double(1 - norm(PR_B(:,1)./(PR_B(:,1) + PR_B(:,2)) - task(2).pT1')^2 / norm(task(2).pT1' - mean(task(2).pT1))^2);
R2TimeB = double(1 - norm( (RT_B(:,1) - task(2).rtc'))^2 / norm(task(2).rtc' - mean(task(2).rtc))^2);
[R2PRU R2PRB R2TimeU R2TimeB]

%%
duration = 250;
x = (1 : length(LIP_HH)) * duration / length(LIP_HH);

figure;
plot(x, LIP_HH','b','DisplayName', 'LIP-Neu','LineWidth',6);
hold on;
plot(x,LIP_HL','-.g','DisplayName', 'LIP-Pre','LineWidth',6);
plot(x,LIP_LH','--r','DisplayName', 'LIP-Nul','LineWidth',6);
 
xlabel('Time (ms)','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
ylabel('Normallized LIP','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
h_legend = legend('show', 'Location', 'SouthEast');
set(h_legend,'FontSize',42);
xlim([0,duration]);
ylim([0,1]);
set(gcf,'paperunits','inches');
set(gcf,'papersize',[10 10]);
set(gcf,'paperposition',[0,0,10, 10]);
set(gca, 'LineWidth',2,'FontWeight','bold');
saveas(gcf, 'LIPResponses_a.fig','fig');
saveas(gcf, 'LIPResponses_a.eps','epsc');

figure;
bias_signal = LIP_HL - LIP_HH;
plot(x,bias_signal,'g');
hold on;
xlabel('Time (ms)','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
ylabel('Difference in LIP','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
xlim([0,duration]);
set(gca, 'LineWidth',2,'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[10 10]);
set(gcf,'paperposition',[0,0,10, 10]);
saveas(gcf, 'LIPResponses_c.fig','fig');
saveas(gcf, 'LIPResponses_c.eps','epsc');

figure;
plot((1: size(LIP_U,2)) * linSqrFit_U(1), mean(LIP_U),'b','LineWidth',6);
hold on;
plot((1: size(LIP_B,2)) * linSqrFit_U(1), mean(LIP_B), '-.g','LineWidth',6) ;
plot((1: size(LIP_B2,2)) * linSqrFit_U(1), mean(LIP_B2), '--r','LineWidth',6) ;
xlabel('Time (ms)','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
ylabel('Normallized LIP','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
ylim([0,max(mean(LIP_B))]);
xlim([0,max(spanU)* linSqrFit_U(1)]);
% set(gca, 'XTick', [-.5 -.25 0 .25 .5],'XMinorTick', 'on', 'LineWidth',2,...
%     'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[10 10]);
set(gcf,'paperposition',[0,0,10, 10]);
saveas(gcf, 'LIPResponses_d.fig','fig');
saveas(gcf, 'LIPResponses_d.eps','epsc');

figure;
%subplot(3,2,4);
%bias_signal = mean(LIP_B-LIP_U);
%bias_signal = LIP_B(1,:) - LIP_U(1,:); 
winSz = 10;
bias_signal = boundary_B - boundary_U;
bias_signal = conv(bias_signal, ones(winSz,1)/winSz,'valid');
bias_signal2 = boundary_B2 - boundary_U;
plot((1:length(bias_signal))*linSqrFit_U(1), bias_signal - min(bias_signal), 'g','LineWidth',6);%, 'DisplayName', 'Model-derived');
hold on;
plot((1:length(bias_signal2))*linSqrFit_U(1), bias_signal2 - min(bias_signal),'--g','LineWidth',6,'DisplayName', 'Neural-derived');
xlabel('Time (ms)','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
ylabel('Dynamic Bias','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
 
xlim([2,length(bias_signal)* linSqrFit_U(1)]);
%ylim([0, 0.25]);
%h_legend = legend('show', 'Location', 'NorthWest');
%set(h_legend,'FontSize',20);
sprintf('Coefficient of Determination R2= \n%.3f (accuracy, neutral),\t %.3f (RT, neutral),\n%.3f (accuracy, bias),\t %.3f (RT, bias)' , R2PRU, R2TimeU, R2PRB, R2TimeB)
set(gcf,'paperunits','inches');
set(gcf,'papersize',[10 10]);
set(gcf,'paperposition',[0,0,10, 10]);
saveas(gcf, 'LIPResponses_e.fig','fig');
saveas(gcf, 'LIPResponses_e.eps','epsc');


figure;
plot(c, PR_U(:,1)./(PR_U(:,1) + PR_U(:,2)), '-','DisplayName', 'Neutral Fit $R^2=0.99$','LineWidth',4);
hold on;
plot(c, PR_B(:,1)./(PR_B(:,1) + PR_B(:,2)), '--g','DisplayName','Bias Fit $R^2=0.98$','LineWidth',4);
plot(c, task(1).pT1,'o',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10, 'DisplayName', 'Neutral Exp');
plot(c, task(2).pT1,'og',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10,'DisplayName', 'Bias Exp');
hold off;
xlim([min(c) max(c)]);
xlabel('Motion Strength','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman')
ylabel('Proportion Correct','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
set(gca, 'XTick', [-.5 -.25 0 .25 .5],'XMinorTick', 'on',  'LineWidth',2,...
    'FontWeight','bold')
h_legend=legend('show','Location','SouthEast');
set(h_legend,'FontSize',24);
set(gcf,'paperunits','inches');
set(gcf,'papersize',[10 10]);
set(gcf,'paperposition',[0,0,10, 10]);
saveas(gcf, 'LIPResponses_g.fig','fig');
saveas(gcf, 'LIPResponses_g.eps','epsc');

figure;
%subplot(3,2,6);
%Show reaction times of correct trials only
plot(c, RT_U(:,1),'-','DisplayName', 'Neutral Fit $R^2=0.97$' ,'LineWidth',4);
hold on;
plot(c, RT_B(:,1),'--g','DisplayName', 'Bias Fit $R^2=0.80$' ,'LineWidth',4);
errorbar(c, task(1).rtc, task(1).rtc_se, 'o', ...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10,'DisplayName', 'Neutral Exp');
errorbar(c, task(2).rtc, task(2).rtc_se, 'go', ...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10,'DisplayName', 'Bias Exp');
hold off;
xlim([min(c) max(c)]);
ylim([350,850]);
xlabel('Motion Strength','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
ylabel('Reaction Time','FontWeight', 'bold','FontSize',42,'FontName','Times New Roman');
h_legend=legend('show','Location','NorthEast');
set(h_legend,'FontSize',24);
set(gca, 'XTick', [-.5 -.25 0 .25 .5],'XMinorTick', 'on', 'LineWidth',2,...
    'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[10 10]);
set(gcf,'paperposition',[0,0,10, 10]);
saveas(gcf, 'LIPResponses_h.fig','fig');
saveas(gcf, 'LIPResponses_h.eps','epsc');
