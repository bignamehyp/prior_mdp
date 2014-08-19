function [psi, t_half] = plotBoundary(dCell,ratio)
colors = 'kbrgmcy';
smoothWinSize = 3;
smoothKernel = ones(1,2 * smoothWinSize + 1) / (2 * smoothWinSize + 1);
figure;
t_half = zeros(length(dCell),1);
subplot(1,2,1);
for i = 1 : length(dCell)
    d = dCell{i};
    t = 1 : size(d,1);
    psi(:,i) = d(:,2) ./ t';
    psi(smoothWinSize + 1:end-smoothWinSize,i) = conv(psi(:,i), smoothKernel, 'valid');
    temp = find(psi(:,i) < 0.75);
    if ~isempty(temp)
        t_half(i) = temp(1);
    else
         t_half(i) = size(psi,1);
    end
    plot(t,psi(:,i), colors(mod(i, length(colors))), 'DisplayName',sprintf('%g',ratio(i)));
    hold on;
end
xlim([1 400]);
ylim([0.5 1]);
xlabel('t', 'FontWeight', 'bold','FontSize',30);
ylabel('\phi^R', 'FontWeight', 'bold','FontSize',30);
legend show;
hold off;
subplot(1,2,2);
plot(ratio,t_half,'-o', 'MarkerSize',10);
xlabel('$\frac{R_P-R_N}{R_S}$', 'FontWeight', 'bold','FontSize',30,'Interpreter','latex');
ylabel('\tau_{1/2}', 'FontWeight', 'bold','FontSize',30);
xlim([100 ratio(end)]);
set(gca, 'XScale', 'log',  'XTick',[10^3 10^6 10^14], 'LineWidth',2,...
    'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[18 12]);
set(gcf,'paperposition',[0,0,18,12]);

saveas(gcf,'boundaries.fig','fig');
saveas(gcf,'boundaries.pdf','pdf');