function [ Prop,  RT, RT_SE, trajectory_R, trajectory_L] = GaussSimulateRT( d, c, mu, sigma, prior)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 5
    prior = 0.5;
end
Prop = zeros(length(c),2);
RT = zeros(length(c),2);
RT_SE = zeros(length(c),2);
trajectory_R = zeros(length(c),size(d,1));
trajectory_L = zeros(length(c),size(d,1));
for i = 1 : length(c)
    [Prop(i,:),  RT(i,:), RT_SE(i,:), trajectory_R(i,:), trajectory_L(i,:)] ...
        = GaussRandomDots(d, c(i), mu, sigma,prior);
end

end

function [Prop, RT, RT_SE, trajectory_R, trajectory_L] =  GaussRandomDots(d, c, mu, sigma, prior)

if nargin < 5
    prior = 0.5;
end

numTrials = 10000;
 
P = 0;
N = 0;
iter = 0;
R_P = 0;
R_N = 0;
STD_P = 0;
STD_N = 0;
nT = size(d,1);
thrR = d(:,2)';
thrL = d(:,1)';

if nargout > 3 
%     trajectories_R = zeros(numTrials,nT);
%     trajectories_L = zeros(numTrials,nT);
    trajectory_R = zeros(1, nT);
    trajectory_L = zeros(1, nT);
end
 
while iter < numTrials
    iter = iter + 1;
    %Generate the observation seequence
    o = mu * c + randn(1,nT) * sigma;
    %Compute the mean observations  
    o = cumsum(o) ./ (1:nT);
    
    tR = find(o >= thrR);
    tL = find(o <= thrL);
    
     
    if isempty(tR)
        tR(1) = nT;
    end

    if isempty(tL)
        tL(1) = nT;
    end
    
    if(tR(1) > 1)
        rB = tR(1) - 1 + (o(tR(1)-1) - thrR(tR(1)-1)) / (thrR(tR(1)) - thrR(tR(1)-1) - o(tR(1)) + o(tR(1)-1) );
    else
        rB = 1;
    end
    if tL(1) > 1
        lB = tL(1) - 1 + (o(tL(1)-1) - thrL(tL(1)-1)) / (thrL(tL(1)) - thrL(tL(1)-1) - o(tL(1)) + o(tL(1)-1));
    else
        lB = 1;
    end
    
    if tR(1) == tL(1)
        if o(nT) < 0;
            tR(1) = tR(1) + 1;
        else
            tL(1) = tL(1) + 1;
        end
    end
    
    if tR(1) < tL(1)
        P = P + 1;
        R_P = R_P + rB;
        STD_P = STD_P + rB * rB;
        if nargout > 3
            tmpR = o;
            tmpR(tR(1):end) = thrR(tR(1):end);
            trajectory_R = trajectory_R + tmpR;
%             tmpR(tR(1):end) = 0;
%             trajectories_R(iter,:) = tmpR;
%             trajectories_L(iter,:) = 0;
        end                 
    end    
     
    
    if tL(1) < tR(1)
        N = N + 1;
        R_N = R_N + lB;
        STD_N = STD_N + lB * lB;
        if nargout > 3
            tmpL = o;
            tmpL(tL(1):end) = thrL(tL(1):end);
            trajectory_L = trajectory_L + tmpL;
%             tmpL(tL(1):end) = 0;
%             trajectories_L(iter,:) = tmpL;
%             trajectories_R(iter,:) = 0;
        end
    end
end
R_P = R_P / P;
R_N = R_N / N;
STD_P = sqrt((STD_P - R_P * R_P * P) / (P - 1));
STD_N = sqrt((STD_N - R_N * R_N * N) / (N - 1));
Prop = [P N];
RT = [R_P R_N];
RT_SE = [STD_P/sqrt(P) STD_N/sqrt(N)];

if nargout > 3    
    trajectory_R = trajectory_R / P;
    trajectory_L = trajectory_L / N;
%    nnz = sum(trajectories_R ~= 0 );
%     cut_off = find(nnz < P/2);
%     nnz(cut_off(1):end) = 0;
%    trajectory_R = sum(trajectories_R) ./ nnz;
%    nnz = sum(trajectories_L ~= 0 );
%     cut_off = find(nnz < N/2);
%     nnz(cut_off(1):end) = 0;
%    trajectory_L = sum(trajectories_L) ./ nnz;
end

end