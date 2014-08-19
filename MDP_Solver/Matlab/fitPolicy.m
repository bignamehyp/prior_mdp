function [Fitted_Curve Est] = fitPolicy(d, deltaT, Bound, isFigure)

if nargin <= 1
    deltaT = 9.5017 * 0.001;
    Bound = 48.6;
    isFigure = false;
end

starting = [8 0.5];
options = '';
phi =  d(:,2);

MAX_PHI = max(phi);
phi = phi / MAX_PHI;
nT = length(phi);
x = (1: nT)';
 
Est = fminsearch(@fitinv, starting, options, x , phi);
Fitted_Curve =  Est(2) * x ./ (x +  Est(1));


if isFigure
figure;
hold off;
plot(deltaT * x, Bound * phi,'ro-');
hold on;
plot(deltaT * x, Bound * Fitted_Curve / MAX_PHI);
xlabel('Time t (s)', 'FontWeight', 'bold','FontSize',30);
ylabel('Policy \pi', 'FontWeight', 'bold','FontSize',30);
hleg = legend('Optimal \phi',sprintf('t_{1/2} = %.3f', deltaT * (Est(1) + sum(phi>=1))),'Location','SouthEast');
set(hleg, 'FontWeight', 'bold','FontSize',30);
set(gca,  'LineWidth',2,...
    'FontWeight','bold');
set(gcf,'paperunits','inches');
set(gcf,'papersize',[12 18]);
set(gcf,'paperposition',[0,0,12,18]);
saveas(gcf, 'decisionBoundary.jpg','jpg');
saveas(gcf, 'decisionBoundary.fig','fig');
end

function sse=fitinv(params,x,y)
tau=params(1);
A = params(2);
Fitted_Curve=   1 - A * x ./(x + tau);
Error_Vector = Fitted_Curve - y;
% When curvefitting, a typical quantity to
% minimize is the sum of squares error
% You could also write sse as
% sse=Error_Vector(:)'*Error_Vector(:);
sse=sum(Error_Vector.^2);

