function [PR RT] = plotPCRT( dCell,ratio)
colors = 'kbrgmcy';
smoothWinSize = 3;
smoothKernel = ones(1,2 * smoothWinSize + 1) / (2 * smoothWinSize + 1);

c = [exp(log(0.01):0.2:log(0.512)) 0.512]';
csize = length(c);
dsize = length(dCell);

PR = zeros(csize,dsize);
RT = zeros(csize,dsize);



for i = 1 : dsize
    d = dCell{i};
    [PR(:,i), RTtemp] = modelRT(d,c);
    RT(:,i) = RTtemp(:,1);
end


figure;
subplot(1,2,1);
for i = 1 : dsize
    plot(c, PR(:,i),colors(i),'DisplayName',sprintf('%g',ratio(i)));
    hold on;
end
hold off;
legend('Location','SouthEast');
xlabel('Coherence','FontWeight', 'bold');
ylabel('Proportion Correct','FontWeight', 'bold');
xlim([0.01, 0.512]);
ylim([0.5 1]);
set(gca, 'XScale', 'log', 'XTick',[0.01 0.1 0.5], 'LineWidth',2,...
    'FontWeight','bold');
subplot(1,2,2);
for i = 1 : dsize
    plot(c, RT(:,i),colors(i),'DisplayName',sprintf('%g',ratio(i)));
    hold on;
end
hold off;
legend show;
xlabel('Coherence', 'FontWeight', 'bold','FontSize',30);
ylabel('POMDP Steps', 'FontWeight', 'bold','FontSize',30);
xlim([0.01, 0.512]);
set(gca, 'XScale', 'log', 'XTick',[0.01 0.1 0.5], 'LineWidth',2,...
    'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[18 12]);
set(gcf,'paperposition',[0,0,18,12]);
saveas(gcf,'PCRTs.fig','fig');
saveas(gcf,'PCRTs.eps','eps');

end

