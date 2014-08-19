function [intercep reci, probit] = reciprobit(RT,linSqrFit_Later, prior)
 %[mu sigma ] = [0.3 0.0455] works for intercep = 6.23

colors = 'kbrgymc';

cmf =  [ 0.5 1 2 5 10 20  50 80 90 95 98 99 99.5  ]*0.01;
%cmf = [1 2 5 10 20 50 80 90 95 98 99] * 0.01;
reci = quantile(RT, cmf); %quantile in descreasing order. 
reci = 1.0 ./ reci;
probit = erfinv(2 * cmf - 1) * sqrt(2);

figure;
starting = [ones(1, length(prior)), 6.2];
for i = 1 : length(prior)
    linFit = [reci(:,i), ones(length(probit),1)] \ probit';
    starting(i) = linFit(1);
end
options = optimset('MaxFunEvals', 1.0e6);
Est = fminsearch(@fitNorm, starting, options, reci, probit');


for i = 2 : length(prior) -1
    plot(reci(:,i), probit , 'LineStyle','none', 'Marker', 'o','Color', colors(mod(i, length(colors)) + 1),  ...
    'DisplayName', sprintf('prior=%.2f', prior(i)));
    %linFit = [reci(:,i), ones(length(probit),1)] \ probit';
    hold on;
    %h2 = plot(reci(:,i), reci(:,i) * linFit(1) + linFit(2) , colors(mod(i, length(colors)) + 1));
    h2 = plot([0; reci(:,i)], [0; reci(:,i)] * Est(i) + Est(end), colors(mod(i, length(colors)) + 1));
    set(get(get(h2,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
end
intercep = Est(end);
legend('show','Location','Best');
hold off;
xlabel('Reciprocal Reaction Time (1/ms)', 'FontWeight','bold', 'FontSize', 30);
%ylim(erfinv(2 * [cmf(1) cmf(end)]  - 1) * sqrt(2));
ylim([erfinv(2 * cmf(1)  - 1) * sqrt(2), Est(end) ]);
ylabel('Probit CDF', 'FontWeight',  'bold', 'FontSize', 30);
annotation(gcf,'textbox',...
    [0.893373923676568 0.0890034360948165 0.0687318489835431 0.0158530760464639],...
    'Interpreter','latex',...
    'String',{'$\infty$'},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');
set(gca, 'XDir', 'reverse', ...
    'LineWidth', 2, 'FontWeight', 'bold', ...
 'XTick',  linSqrFit_Later(1)./([500 300 200 180] - linSqrFit_Later(2)), 'XTickLabel', [ 500 300 200 180], ...
 'YTick', erfinv(2 *[0.001 0.01 0.1 0.5 0.9 0.99 0.999] -1 ) * sqrt(2), ...
  'YTickLabel', [0.001 0.01 0.1 0.5 0.9 0.99 0.999]);
end


function sse=fitNorm(params,x,y)
sse = 0;

for i = 1 : length(params) - 1
    sse = sse + norm(params(i) * x(:,i) + params(end) - y)^2;
end

end

