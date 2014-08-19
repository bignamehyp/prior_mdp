function  plotPolicy( d , isSave)
if nargin <= 1
    isSave = 0;
end

maxBoundary = 3;



%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nT = size(d,1);
yResolution = 1.0e5;

 

policy = zeros(yResolution,nT);

leftB = zeros(1,nT);
rightB = leftB;

if d(end,2) < 1 %GaussPolicy format
    leftB = d(:,1)/maxBoundary/2 * yResolution + yResolution/2;
    rightB = d(:,2)/maxBoundary/2 * yResolution + yResolution/2;
else %BetaPolicy format
    leftB = d(:,1) ./(1:nT)' * yResolution;
    rightB = d(:,2) ./(1:nT)' * yResolution;
end

leftB(leftB < 1) = 1;
leftB(leftB > yResolution) = yResolution;
rightB(rightB < 1) = 1;
rightB(rightB > yResolution) = yResolution;

% 
% for t = 1 : nT
%     leftB(t) = (1 -  betacdf(0.5, d(t,1) +  0.01, t - d(t,1) - 0.01)) ...
%             * yResolution;
% 
%     rightB(t) = (1 -  betacdf(0.5, d(t,2) - 0.01, t - d(t,2) + 0.01)) ...
%             * yResolution;    
% end


window = 2;
leftB(window:end) = conv(leftB, ones(window,1)/window,'valid');
rightB(window:end) = conv(rightB, ones(window,1)/window,'valid');
leftB = floor(leftB)+1;
rightB = floor(rightB);

for t  = 1: nT
    policy(1 : (leftB(t) - 1),t) = 1; %left
    policy(leftB(t) : rightB(t),t) = 3 ; %Sample
    policy((rightB(t)+1):end,t) = 2; %right;
end


figure;
image( 1:nT, (1:yResolution)/yResolution, policy * 16 );
ylabel('$\bar{s}$','FontWeight', 'bold','FontSize',42,'Interpreter','latex');
xlabel('Time Step t','FontWeight', 'bold','FontSize',30);
set(gca, 'LineWidth',2,'FontSize', 30,'YDir','normal', ...
    'FontWeight','bold', 'YTick',[0.25 0.5 0.75],'YTickLabel',[-1.5 0 1.5]);
set(gca, 'Units', 'normalized','Position', [0.16 0.18 0.8 0.8]);
set(gcf,'paperunits','inches');
set(gcf,'papersize',[12 8]);
set(gcf,'paperposition',[0,0,12,8]);
if isSave
 saveas(gcf,'GaussPolicy.fig','fig');
 saveas(gcf,'GaussPolicy.pdf','pdf');
end
end

