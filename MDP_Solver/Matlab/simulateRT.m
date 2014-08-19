function [ Prop,  RT, RT_STD, LIP] = simulateRT( d, c, b0, nT, w_prior, option)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Prop = zeros(length(c),2);
RT = zeros(length(c),2);
RT_STD = zeros(length(c),2);
LIP = zeros(length(c),nT);

for i = 1 : length(c)
    [Prop(i,:),  RT(i,:), LIP(i,:) RT_STD(i,:)] = randomDots(b0,c(i),nT,d,w_prior);
end

if(option > 0)
    %Cover the probability to coherence.
    
    %Fit the proportion correct / reaction time to the DDM model
    %where A = Estimates(1) and k = Estimates(2)
    %RT(c) = A/(kc) tanh(A*k*c)
    %We estimate A and c from the above equation
    %With the same set of paramters,
    %Predict the proportion correct
    %PC(c) = 1 / ( 1 + exp(-2Ak|c|) );
    
    starting = rand(1,2);
    options = '';
    Estimates1 = fminsearch(@fitDDM, starting, options, c, RT(:,1));
    Fitted_Curve= Estimates1(1)/Estimates1(2)./c .* tanh(Estimates1(1)*Estimates1(2)*c);
    PropCorrect1 = 1.0 ./ ( 1 + exp( - 2 * abs(Estimates1(1) * Estimates1(2) * c)));
    
    % starting = rand(1,2);
    % options = optimset('Display','iter');
    % Estimates2 = fminsearch(@fitDDM, starting, options, c, RT(:,2));
    % Fitted_Curve2 = Estimates2(1)/Estimates2(2)./c .* tanh(Estimates2(1)*Estimates2(2)*c);
    % PropCorrect2 = 1.0 ./ ( 1 + exp( - 2 * abs(Estimates2(1) * Estimates2(2) * c)));
    
    figure;
    subplot(2,1,1);
    errorbar(c, RT(:,1), RT_STD(:,1)./sqrt(Prop(:,1)),'ob');
    hold on;
    errorbar(c, RT(:,2), RT_STD(:,2)./sqrt(Prop(:,2)),'og');
    plot(c, Fitted_Curve,'-b');
    %plot(c, Fitted_Curve2,'-g');
    legend('Reaction time for correct choices \pm SEM', 'Reaction time for incorrect choices \pm SEM',...
        'Location','SouthWest');
    xlabel('Coherence');
    ylabel('Reaction Time');
    xlim([c(1), c(end)]);
    set(gca, 'XScale','log','XTick',[0 0.032 0.064 0.128 0.256 0.512]);
    subplot(2,1,2);
    semilogx(c, Prop(:,1) ./ (Prop(:,1) + Prop(:,2)),'o');
    hold on;
    semilogx(c, PropCorrect1,'-');
    xlabel('Coherence');
    ylabel('Proportion Correct');
    xlim([c(1), c(end)]);
    set(gca, 'XTick',[0 0.032 0.064 0.128 0.256 0.512]);
end
if(option > 1)
    set(gcf,'paperunits','inches');
    set(gcf,'papersize',[11 11]);
    set(gcf,'paperposition',[0,0,11,11]);
    saveas(gcf,'PCRT.fig','fig');
    saveas(gcf,'PCRT.jpg','jpg');
end
end

