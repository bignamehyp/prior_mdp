%% Plotting setting
%clear;
set(0, 'DefaultAxesFontSize',30);
set(0, 'DefaultLineLineWidth',6);


%% Speed and Accuracy Trade-off Palmer

c = [0.008, 0.016,0.032,0.064,0.128,0.256,0.512];

PR0_JP = [.55 .62 .68 .77 .94 .99  1 ]';
PR0_JP_SE = [0.05 0.04 0.04 0.04 0.02 0.01 0]';
RT0_JP = [502 515 494 480 457 372 356]';
RT0_JP_SE = [14 14 13 12 9 5 3 ]';

PR1_JP = [.61 .68 .82 .96 .98 1 1]';
PR1_JP_SE = [0.05 0.04 0.04 0.02 0.01 0 0]';
RT1_JP = [1016 1011 981 853 627 465 402]';
RT1_JP_SE = [38 42 32 30 19 8 4 ]';

PR2_JP = [0.63 0.6 0.81 0.97 0.99 1 1  ]';
PR2_JP_SE = [0.05 0.05 0.04 0.01 0.01 0 0]';
RT2_JP = [1890 1927 1530 1220 710 487 412]';
RT2_JP_SE = [106 101 72 58 24 10 5]';

R_P_Set = [1000:100:2000];
mu_Set = 1.5;
sigma_Set = 1;

%[ R_P_min, mu_min, sigma_min, energy_min ] = fitPRRT(PR1_JP, RT1_JP, c, R_P_Set, mu_Set, sigma_Set)


d0 = load('Policy/GaussPolicy_-0.1_60.0_0.0.txt');
d1 = load('Policy/GaussPolicy_-0.1_1500.0_0.0.txt');
d2 = load('Policy/GaussPolicy_-0.1_1000.0_0.0.txt');


nT = round(size(d0, 1)* 0.8);

[PR0, RT0] = GaussSimulateRT(d0(1:nT,:),c,1.5,1);
[PR1, RT1] = GaussSimulateRT(d1(1:nT,:),c,1.5,1);
[PR2, RT2] = GaussSimulateRT(d2(1:nT,:),c,1.5,1);

RT0(c < 0,:) = fliplr(RT0(c < 0,:));
RT1(c < 0,:) = fliplr(RT1(c < 0,:));
RT2(c < 0,:) = fliplr(RT2(c < 0,:));


% linSqrFit =  [[RT0(:,1);RT1(:,1); RT2(:,2)], ones(length(c) * 3,1)] ...
%     \ [RT0_JP; RT1_JP; RT2_JP];

linSqrFit0 =  [RT0(:,1), ones(length(c),1)] \ RT0_JP;
linSqrFit1 =  [RT1(:,1), ones(length(c),1)] \ RT1_JP;
linSqrFit2 =  [RT2(:,1), ones(length(c),1)] \ RT2_JP;


RT0 = RT0 * linSqrFit0(1) + linSqrFit0(2);

RT1 = RT1 * linSqrFit1(1) + linSqrFit1(2);

RT2 = RT2 * linSqrFit2(1) + linSqrFit2(2);

 
figure;
subplot(1,2,1);
plot(c, PR0(:,1) ./ (PR0(:,1) + PR0(:,2)));
hold on;
plot(c, PR1(:,1) ./ (PR1(:,1) + PR1(:,2)),'--r');
% plot(c, PR2(:,1) ./ (PR2(:,1) + PR2(:,2)),'--g');
errorbar(c, PR0_JP, PR0_JP_SE, 'xb',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',25);
errorbar(c, PR1_JP, PR1_JP_SE, 'xr',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',25);
% errorbar(c, PR2_JP, PR2_JP_SE, 'xg',...
%     'MarkerFaceColor',[0 0 0],'MarkerSize',25);
xlabel('Coherence','FontWeight', 'bold');
ylabel('Proportion Correct','FontWeight', 'bold');
xlim([0 0.512]); 
ylim([0.5 1]);
set(gca, 'XScale', 'log','LineWidth',2,'FontWeight','bold');

subplot(1,2,2);
plot(c, RT0(:,1));
hold on;
plot(c, RT1(:,1),'--r');
% plot(c, RT2(:,1),'--g');
errorbar(c,RT0_JP, RT0_JP_SE,'xb',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',25);
errorbar(c, RT1_JP, RT1_JP_SE, 'xr',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',25);
% errorbar(c, RT2_JP, RT2_JP_SE, 'xg',...
%     'MarkerFaceColor',[0 0 0],'MarkerSize',25);
xlabel('Coherence', 'FontWeight', 'bold','FontSize',30);
ylabel('Reaction Time (ms)', 'FontWeight', 'bold','FontSize',30);
xlim([0, 0.512]); 
set(gca, 'xScale', 'log','xScale', 'log','YScale', 'log', 'LineWidth',2,...
    'FontWeight','bold');


%% Effects of Prior Probability on reaction times


% Coefficient of Determination
% Subject SK, speed
% 0.997 (accuracy, neutral), 0.943 (RT, neutral),
% 0.989 (accuracy, bias),   0.942 (RT, bias)
% Subject LH, speed
% 0.988 (accuracy, neutral), 0.955 (RT, neutral),
% 0.984 (accuracy, bias),     0.886 (RT, bias)
% Monkey (prior = 0.8)
% 0.999 (accuracy, neutral), 0.987 (RT, neutral),
% 0.979 (accuracy, bias), 0.846 (RT, bias)
% Monkey (prior = 0.7)
% 0.999 (accuracy, neutral), 0.987 (RT, neutral),
% 0.980 (accuracy, bias),     0.977 (RT, bias)


%[mu sigma] = (1.45 1.1) Ssp
%[mu sigma] = (1.4 1.05) monk prior = 0.7
%[mu sigma] = (1.45 0.95) Lsp

% d0_B = load('Policy/GaussPolicy_biased_0.7.txt');
% d0_U = load('Policy/GaussPolicy_-0.1_100.0_0.0.txt');
 
d0_B = load('Policy/GaussPolicy_0.70_-0.1_100.0_0.0_200.txt'); %mu sigma = 1 1 
d0_U = load('Policy/GaussPolicy_0.50_-0.1_100.0_0.0_200.txt');


nT = size(d0_B,1);

boundary_U = stateToProb(d0_U(1:nT,2)', sigma./sqrt(1:nT), 0.5);
boundary_U = log10(boundary_U ./(1 - boundary_U));
boundary_B =  boundary_U - 0.4 *(1 - exp(-(1:nT)/nT*16)); 
cdf = 1.0 ./ (1.0 + 10.^(boundary_B- log10(prior / ( 1 - prior) ))); 
cdf = erfinv(1 - 2 * cdf) * sqrt(2) * sigma./sqrt(1:nT);
d0_B(1:nT,1) = d0_B(1:nT,1) + d0_B(1:nT,2) - cdf';
d0_B(1:nT,2) = cdf;

 

filename = 'monk';
load(sprintf('%s_data.mat',filename));
c = task(1).coh_set;

%%%[mu, sigma, Energy] = fitBiasData(d0_U, task);

mu = 1;
sigma = 1;


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

sprintf('Coefficient of Determination R2= \n%.3f (accuracy, neutral),\t %.3f (RT, neutral),\n%.3f (accuracy, bias),\t %.3f (RT, bias)' , R2PRU, R2TimeU, R2PRB, R2TimeB)
figure;
subplot(2,1,1);
plot(c, PR_U(:,1)./(PR_U(:,1) + PR_U(:,2)), '-','DisplayName', 'Neutral Fit');
hold on;
plot(c, PR_B(:,1)./(PR_B(:,1) + PR_B(:,2)), '--g','DisplayName','Bias Fit');
plot(c, task(1).pT1,'o',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10, 'DisplayName', 'Neutral Exp');
plot(c, task(2).pT1,'og',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10,'DisplayName', 'Bias Exp');
hold off;
xlim([min(c) max(c)]);
ylabel('Proportional Correct','FontWeight', 'bold','FontSize',30);
set(gca, 'XTick', [-.5 -.25 0 .25 .5],'XMinorTick', 'on',  'LineWidth',2,...
    'FontWeight','bold')
%legend('show','Location','SouthEast');

subplot(2,1,2);
%Show reaction times of correct trials only
plot(c, RT_U(:,1),'-','DisplayName', 'Neutral Fit');
hold on;
plot(c, RT_B(:,1),'--g','DisplayName', 'Bias Fit');
errorbar(c, task(1).rtc, task(1).rtc_se, 'o', ...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10,'DisplayName', 'Neutral Exp');
errorbar(c, task(2).rtc, task(2).rtc_se, 'go', ...
    'MarkerFaceColor',[0 0 0],'MarkerSize',10,'DisplayName', 'Bias Exp');
hold off;
xlim([min(c) max(c)]);
xlabel('Motion Strength','FontWeight', 'bold','FontSize',30);
ylabel('Reaction Time','FontWeight', 'bold','FontSize',30);
%legend('show');
set(gca, 'XTick', [-.5 -.25 0 .25 .5],'XMinorTick', 'on', 'LineWidth',2,...
    'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[9 15]);
set(gcf,'paperposition',[0,0,9,15]);
% saveas(gcf,sprintf('GaussPrior_%s.fig',filename),'fig');
% saveas(gcf,sprintf('GaussPrior_%s.jpg',filename),'jpg');

%% Fixed Time Version Prior
c_s = [-0.32 -0.16 -0.08 -0.04 0 0.04 0.08 0.16 0.32];
PR_Nul_Short = [0.0333 0.0929 0.1365 0.1767 0.2287 0.3508 0.4391 0.5463 0.7943];
PR_Neu_Short = [0.0974 0.2233 0.3195 0.4109 0.4572 0.5530 0.6248 0.7875 0.9321];
PR_Pre_Short = [0.2909 0.4560 0.5878 0.7032 0.7240 0.8030 0.8547 0.9038 0.9708];
PR_Nul_Short_SE = [0.0613 0.0965 0.1260 0.1282 0.1312 0.2083 0.2621 0.2791 0.2373] / sqrt(72);
PR_Neu_Short_SE = [0.098 0.1471 0.1759 0.1703 0.1427 0.1770 0.1589 0.1414 0.08] / sqrt(72);
PR_Pre_Short_SE = [0.2317 0.2463 0.2435 0.1806 0.1493 0.1167 0.1247 0.1028 0.0510] / sqrt(72);



c_l = [-0.24 -0.12 -0.06 -0.03 0 0.03 0.06 0.12 0.24];
PR_Nul_Long = [0.0081 0.0842 0.2063 0.2471 0.3319 0.4447 0.5970 0.8004 0.9720];
PR_Neu_Long = [0.019 0.1489 0.2644 0.3555 0.4120 0.5267 0.6801 0.8698 0.9840];
PR_Pre_Long = [0.117 0.2923 0.4753 0.5011 0.6428 0.7322 0.8441 0.9535 0.9968];
PR_Nul_Long_SE = [0.0263 0.1041 0.1338 0.1490 0.1630 0.2091 0.1978 0.1856 0.0750] / sqrt(44);
PR_Neu_Long_SE = [0.0351 0.1153 0.1752 0.1578 0.1331 0.1559 0.1542 0.1033 0.0369] / sqrt(44);
PR_Pre_Long_SE = [0.1565 0.2214 0.1912 0.1806 0.1804 0.1763 0.1309 0.0631 0.0148] / sqrt(44);

%linSqrFit = [c_s', ones(length(c_s),1) ] \ (sqrt(2) * erfinv( 2 * PR_Pre_Short - 1)');

mu = 1;
sigma = 1;

% R_P_H = normcdf(linSqrFit(2));
% R_P_L = 1 - R_P_H;
% R_N = 0;

prior = 0.80;
d0_B = load('Policy/GaussPolicy_0.80_-0.1_100.0_0.0_30.txt');
d0_U = load('Policy/GaussPolicy_0.50_-0.1_100.0_0.0_30.txt');

PR = GaussSimulateRT( d0_U, c_s, mu, sigma);
PR_HH_Short = PR(:,1) ./ (PR(:,1) + PR(:,2));
PR = GaussSimulateRT( d0_B, c_s, mu, sigma, prior);
PR_HL_Short = PR(:,1) ./ (PR(:,1) + PR(:,2));
PR_LH_Short = flipud(1 - PR_HL_Short);

% PR_HH_Short = fixedDurationRandomDots(c_s, mu, sigma, R_P_H, R_P_H, R_N);
% PR_HL_Short = fixedDurationRandomDots(c_s, mu, sigma, R_P_H, R_P_L, R_N);
% PR_LH_Short = fixedDurationRandomDots(c_s, mu, sigma, R_P_L, R_P_H, R_N);



R2PRHH = double(1 - norm(PR_HH_Short'  - PR_Neu_Short)^2 / norm(PR_Neu_Short - mean(PR_Neu_Short))^2);
R2PRHL = double(1 - norm(PR_HL_Short'  - PR_Pre_Short)^2 / norm(PR_Pre_Short - mean(PR_Pre_Short))^2);
R2PRLH = double(1 - norm(PR_LH_Short'  - PR_Nul_Short)^2 / norm(PR_Nul_Short - mean(PR_Nul_Short))^2);

[R2PRHH R2PRHL R2PRLH]

figure;
plot(c_s, PR_HH_Short, '-b','DisplayName', 'Neu');
hold on;
plot(c_s, PR_HL_Short, '-k','DisplayName', 'Pref');
plot(c_s, PR_LH_Short, '-g','DisplayName', 'Null');
errorbar(c_s, PR_Neu_Short, PR_Neu_Short_SE, 'ob','DisplayName', 'Neu');
errorbar(c_s, PR_Pre_Short, PR_Pre_Short_SE, 'ok','DisplayName', 'Pref');
errorbar(c_s, PR_Nul_Short, PR_Nul_Short_SE, 'og','DisplayName', 'Null');
xlim([min(c_s) max(c_s)]);
xlabel('Motion Strength');
ylabel('Proportional Correct');
legend_handle = legend('show', 'Location', 'NorthWest');
set(legend_handle, 'Box', 'off')


%linSqrFit = [c_l', ones(length(c_l),1) ] \ (sqrt(2) * erfinv( 2 * PR_Pre_Long - 1)');
%mu = linSqrFit(1);
figure
mu = 1;
sigma = 1;
prior = 0.66;
d0_B = load('Policy/GaussPolicy_0.66_-0.1_100.0_0.0_100.txt');
d0_U = load('Policy/GaussPolicy_0.50_-0.1_100.0_0.0_100.txt');
PR = GaussSimulateRT( d0_U, c_l, mu, sigma);
PR_HH_Long = PR(:,1) ./ (PR(:,1) + PR(:,2));
PR = GaussSimulateRT( d0_B, c_l, mu, sigma, prior);
PR_HL_Long = PR(:,1) ./ (PR(:,1) + PR(:,2));
PR_LH_Long = flipud(1 - PR_HL_Long);


R2PRHH = double(1 - norm(PR_HH_Long'  - PR_Neu_Long)^2 / norm(PR_Neu_Short - mean(PR_Neu_Short))^2);
R2PRHL = double(1 - norm(PR_HL_Long'  - PR_Pre_Long)^2 / norm(PR_Pre_Short - mean(PR_Pre_Short))^2);
R2PRLH = double(1 - norm(PR_LH_Long'  - PR_Nul_Long)^2 / norm(PR_Nul_Short - mean(PR_Nul_Short))^2);
[R2PRHH R2PRHL R2PRLH]

% R_P_H = normcdf(linSqrFit(2));
% R_P_L = 1 - R_P_H;
% R_N = 0;
% 
% PR_HH_Long = fixedDurationRandomDots(c_s, mu, sigma, R_P_H, R_P_H, R_N);
% PR_HL_Long = fixedDurationRandomDots(c_s, mu, sigma, R_P_H, R_P_L, R_N);
% PR_LH_Long = fixedDurationRandomDots(c_s, mu, sigma, R_P_L, R_P_H, R_N);

figure;
plot(c_l, PR_HH_Long, '-b','DisplayName', 'Neu');
hold on;
plot(c_l, PR_HL_Long, '-k','DisplayName', 'Pref');
plot(c_l, PR_LH_Long, '-g','DisplayName', 'Null');
errorbar(c_l, PR_Neu_Long, PR_Neu_Long_SE, 'ob','DisplayName', 'Neu');
errorbar(c_l, PR_Pre_Long, PR_Pre_Long_SE, 'ok','DisplayName', 'Pref');
errorbar(c_l, PR_Nul_Long, PR_Nul_Long_SE, 'og','DisplayName', 'Null');
xlim([min(c_l) max(c_l)]);
xlabel('Motion Strength');
ylabel('Proportional Correct');
legend_handle = legend('show', 'Location', 'NorthWest');
set(legend_handle, 'Box', 'off')

set(gcf,'papersize',[12 6]);
set(gcf,'paperposition',[0,0,12,6]);
saveas(gcf,'priorFixDuration.fig','fig');
 


%% 

%Assymetric Rewards  Fixed Duration 
% |mu*c/sigma| < 8
c = [-0.48 -0.24 -0.12 -0.06 -0.03 -0.015 0 0.015, 0.03, 0.06, 0.12 0.24 0.48];

PR_LL_Monkey = 0.01 * [0.4412    5.7663   20.6631   34.8685   45.3458   48.1665   53.9602   57.6485   63.2149   67.7702   82.3921   96.1181   99.9265];
PR_LH_Monkey = 0.01 * [0.0735    0.5277    3.4024    8.2395   10.7901   10.8187   14.1848   17.8686   19.5175   27.0115   45.4243   76.6338   97.3102];
PR_HL_Monkey = 0.01 * [2.4496   25.7552   62.1102   75.5467   80.9535   86.6039   86.6190   90.4944   93.2612   94.7494   96.9185   99.4363   99.9999];
PR_HH_Monkey = 0.01 * [0.2206    4.6174   19.5695   33.7645   42.7872   44.1819   51.2724   54.7453   57.9309   67.6536   82.2595   96.0455   99.9265];

PR_LL_STD = 0.01 * [0.1966    1.0943    1.6778    1.8092    1.6140    2.1540    2.1208    2.2010    2.0677    2.0939    1.8113    0.9258    0.0735];
PR_LH_STD = 0.01 * [0.0735    0.2094    0.6219    0.9549    1.1360    1.1878    1.2751    1.6751    1.8419    2.3344    2.8657    3.1222    0.9421];
PR_HL_STD = 0.01 * [0.6095    2.4289    2.3593    1.9986    1.7044    1.3563    1.6757    1.2895    0.9652    0.8641    1.0018    0.2696    0.0];
PR_HH_STD = 0.01 * [0.2206    0.9020    1.9081    2.1115    2.0891    2.2881    2.4434    1.9795    2.5838    2.2621    1.9569    1.0150    0.0735];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PR = normcdf( (c * mu + theta), 0, 1);
% c * mu + theta = probit(PR) = sqrt(2) * erfinv(2 * PR - 1);
% theta = probit(R_P_H/ (R_P_H + R_P_L))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% linSqrFit = [c', ones(length(c),1) ] \ (sqrt(2) * erfinv( 2 * PR_HL_Monkey - 1)');
% 
% mu = linSqrFit(1);
mu = 1.4;
sigma = 1;

% 
% R_P_H = normcdf(linSqrFit(2));
% R_P_L = 1 - R_P_H;
% R_N = 0;
% 
% PR_HH = fixedDurationRandomDots(c, mu, sigma, R_P_H, R_P_H, R_N);
% PR_HL = fixedDurationRandomDots(c, mu, sigma, R_P_H, R_P_L, R_N);
% PR_LH = fixedDurationRandomDots(c, mu, sigma, R_P_L, R_P_H, R_N);
% PR_LL = fixedDurationRandomDots(c, mu, sigma, R_P_L, R_P_L, R_N);
prior = 0.9;
d0_B = load('Policy/GaussPolicy_0.90_-0.1_100.0_0.0_30.txt');
d0_U = load('Policy/GaussPolicy_0.50_-0.1_100.0_0.0_30.txt');
PR = GaussSimulateRT( d0_U, c, mu, sigma);
PR_HH = PR(:,1) ./ (PR(:,1) + PR(:,2));

d0_U = load('Policy/GaussPolicy_0.50_-0.1_50.0_0.0_30.txt');
PR = GaussSimulateRT( d0_U, c, mu, sigma);
PR_LL = PR(:,1) ./ (PR(:,1) + PR(:,2));

PR = GaussSimulateRT( d0_B, c, mu, sigma, prior);
PR_HL = PR(:,1) ./ (PR(:,1) + PR(:,2));
PR_LH = flipud(1 - PR_HL);




R2PRHH = double(1 - norm(PR_HH'  - PR_HH_Monkey)^2 / norm(PR_HH_Monkey - mean(PR_HH_Monkey))^2);
R2PRLL = double(1 - norm(PR_LL'  - PR_LL_Monkey)^2 / norm(PR_LL_Monkey - mean(PR_LL_Monkey))^2);

R2PRHL = double(1 - norm(PR_HL'  - PR_HL_Monkey)^2 / norm(PR_HL_Monkey - mean(PR_HL_Monkey))^2);
R2PRLH = double(1 - norm(PR_LH' - PR_LH_Monkey)^2 / norm(PR_LH_Monkey - mean(PR_LH_Monkey))^2);

[R2PRHH R2PRLL R2PRHL R2PRLH]


figure;
 plot(c, PR_HL, '-k','DisplayName', 'HL-model');
hold on;
plot(c, PR_LH, '-g','DisplayName', 'LH-model');
plot(c, PR_HH, '-r','DisplayName', 'HH-model');
plot(c, PR_LL, '-b','DisplayName', 'LL-model');
errorbar(c, PR_HL_Monkey, PR_HL_STD, 'ok','DisplayName', 'HL-monkey');
errorbar(c, PR_LH_Monkey, PR_LH_STD, 'og','DisplayName', 'LH-monkey');
errorbar(c, PR_HH_Monkey, PR_HH_STD, 'or','DisplayName', 'HH-monkey');
errorbar(c, PR_LL_Monkey, PR_LL_STD, 'ob','DisplayName', 'LL-monkey');


xlim([min(c) max(c)]);
xlabel('Motion Strength','FontWeight', 'bold','FontSize',30);
ylabel('Proportional Correct','FontWeight', 'bold','FontSize',30);
legend_handle = legend('show', 'Location', 'NorthWest');
set(legend_handle, 'Box', 'off')
 
set(gca,  'LineWidth',2,'FontWeight','bold')
set(gcf,'papersize',[12 9]);
set(gcf,'paperposition',[0,0,12,9]);
saveas(gcf,'fixDuration.fig','fig');
 
% Assymetric Rewards reaction time version
% sigma = 1;
% dHL = load('GaussPolicy_-0.1_100.0_50.0_0.0_0.0.txt');
% dLH = load('GaussPolicy_-0.1_50.0_100.0_0.0_0.0.txt');
% dHH = load('GaussPolicy_-0.1_100.0_0.0.txt');
% dLL = load('GaussPolicy_-0.1_50.0_0.0.txt');
% 
% 
% 
% [PR_HL, RT_HL] = GaussSimulateRT(dHL,c,mu, sigma);
% [PR_LH, RT_LH] = GaussSimulateRT(dLH,c,mu, sigma);
% [PR_HH, RT_HH] = GaussSimulateRT(dHH,c,mu, sigma);
% [PR_LL, RT_LL] = GaussSimulateRT(dLL,c,mu, sigma);
% 
% RT_HL(c < 0,:) = fliplr(RT_HL(c < 0,:));
% RT_LH(c < 0,:) = fliplr(RT_LH(c < 0,:));
% RT_HH(c < 0,:) = fliplr(RT_HH(c < 0,:));
% RT_LL(c < 0,:) = fliplr(RT_LL(c < 0,:));
% 
% figure;
% subplot(2,1,1);
% plot(c, PR_HL(:,1)./(PR_HL(:,1) + PR_HL(:,2)), '-','DisplayName', 'HL');
% hold on;
% plot(c, PR_LH(:,1)./(PR_LH(:,1) + PR_LH(:,2)), '--g','DisplayName','LH');
% plot(c, PR_HH(:,1)./(PR_HH(:,1) + PR_HH(:,2)), '--r','DisplayName','HH');
% plot(c, PR_LL(:,1)./(PR_LL(:,1) + PR_LL(:,2)), '--c','DisplayName','LL');
% hold off;
% xlim([min(c) max(c)]);
% ylabel('Proportional Correct','FontWeight', 'bold','FontSize',30);
% legend('show');
% subplot(2,1,2);
% plot(c, RT_HL(:,1),'-','DisplayName', 'HL');
% hold on;
% plot(c, RT_LH(:,1),'--g','DisplayName', 'LH');
% plot(c, RT_HH(:,1),'--r','DisplayName', 'HH');
% plot(c, RT_LL(:,1),'--c','DisplayName', 'LL');
% xlim([min(c) max(c)]);
% legend('show');
% xlabel('Motion Strength','FontWeight', 'bold','FontSize',30);
% ylabel('Reaction Time','FontWeight', 'bold','FontSize',30);
% set(gca,  'LineWidth',2,'FontWeight','bold')

%%
%Later Model;
d5 = load('Policy/GaussPolicy_-0.1_100.0_0.0.txt');
d6 = load('Policy/GaussPolicy_biased_0.75.txt');
d7 = load('Policy/GaussPolicy_biased_0.9.txt');
d8 = load('Policy/GaussPolicy_biased_0.95.txt');

 %[mu sigma ] = [0.3 0.0457] works for intercep = 6.20
sigma = 0.0457;
mu = 0.3;
nTrials = 10000;
% nT = size(d5,1);
% o = mu * ones(1,nT);
% o(1) = 0;
% o = cumsum(o) ./ (1:nT);
% o(o' > d8(:,2)) = [];
% plotPolicy(d8);
% hold on;
% plot(o/6 + 0.5,'k');
% hold off;
% ylim([0.33,0.67]);
% set(gca, 'LineWidth',2,'FontSize',48,...
%      'FontWeight','bold', 'YTick',[0.33 0.5 0.67],'YTickLabel',[-1 0 1]);
% saveas(gcf,'GaussPolicy.fig','fig');
% saveas(gcf,'GaussPolicy.jpg','jpg');

prior = [0.05 0.1 0.25 0.5, 0.75, 0.9, 0.95];

LaterRT = zeros(nTrials,length(prior));

LaterRT(:,1) = GaussLaterRT(d8,sigma,-mu,nTrials);
LaterRT(:,2) = GaussLaterRT(d7,sigma,-mu,nTrials);
LaterRT(:,3) = GaussLaterRT(d6,sigma,-mu,nTrials);
LaterRT(:,4) = GaussLaterRT(d5,sigma,mu,nTrials);
LaterRT(:,5) = GaussLaterRT(d6,sigma,mu,nTrials);
LaterRT(:,6) = GaussLaterRT(d7,sigma,mu,nTrials);
LaterRT(:,7) = GaussLaterRT(d8,sigma,mu,nTrials);


medianRT = median(LaterRT);
linSqrFit_Later =  [medianRT', ones(length(prior),1)] \ [273;250;241;202;187;179;176];

realLaterRT = LaterRT*linSqrFit_Later(1) + linSqrFit_Later(2);

intercep = reciprobit(LaterRT, linSqrFit_Later, prior)
set(gcf,'paperunits','inches');
set(gcf,'papersize',[15 8]);
set(gcf,'paperposition',[0,0,15,8]);
% saveas(gcf,'GaussReciprobit.fig','fig');
% saveas(gcf,'GaussReciprobit','jpg');

% plot the median RT vs prior figure   
medianRealRT = medianRT*linSqrFit_Later(1) + linSqrFit_Later(2);
linMedian = [log(prior'), ones(length(prior),1)]\ medianRealRT';
linExpMedian = [log(prior'), ones(length(prior),1)] \ [273;250;241;202;187;179;176];

figure;
semilogx(prior, medianRealRT,'*','MarkerFaceColor',[0 0 0],'MarkerSize',12);
hold on;
h2=semilogx(prior, log(prior) * linMedian(1) + linMedian(2));
set(get(get(h2,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
semilogx(prior, [273;250;241;202;187;179;176], '*k',...
    'MarkerFaceColor',[0 0 0],'MarkerSize',12);
h4 = semilogx(prior, ...
    log(prior) * linExpMedian(1) + linExpMedian(2),'--k');
set(get(get(h4,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
legend(sprintf('Model Prediction \n y = %.2f x + %.2f', linMedian(1), linMedian(2)), ...
    sprintf('Exp. Data \n y = %.2f x + %.2f', linExpMedian(1), linExpMedian(2)));
xlim([0.04 0.96]);
xlabel('Log Prior','FontWeight', 'bold','FontSize',30);
ylabel('Median RT(ms)','FontWeight', 'bold','FontSize',30);
set(gca, 'XTick', [0.1 0.2 0.5 0.9],'XMinorTick', 'on', 'LineWidth',2,...
    'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[15 8]);
set(gcf,'paperposition',[0,0,15,8]);
saveas(gcf,'GaussLaterMedian.fig','fig');
saveas(gcf,'GaussLaterMedian','jpg');
