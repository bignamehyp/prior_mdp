function [F, L] = decisionPOMDP(c, r_p, r_n, r_s)
%Inputs:
%       c, motion strength
%       r_p, reward for making the correct action
%       r_n, reward for making the incorrect action
%       r_s, reward for sampling
if nargin < 1
    c = 0.3;
end
if nargin < 2
    r_p = 20;
end
if nargin < 3
    r_n = -100;
end
if nargin < 4
    r_s = -1;
end

nb = 11;

nu = 602;

F = ones(nb,nb, nu); %Transition probability of belief state
L = zeros(nb,nu); %Reward function R(b,a) = \sum_s b(s) r(s,a)
T = zeros(nb,nb);
b =  linspace(0,1,nb);
OR = (c+1)/2.0; %Likelihood function = P(o=O_R | s=S_R) 
OL = 1 - OR;
for i = 1 : nb %P(s = S_R)   
    %Probability that the next belief is b2R
    P_b2R = OR * b(i) + OL * (1 - b(i));
    %Probability that the next belief is b2L
    P_b2L = OL * b(i) + OR * (1 - b(i));
    %Action = A_S sampling
    if P_b2R ~= 0
        b2R = OR * b(i)/P_b2R; %The next possible belief when o_t+1 = O_R
    else
        b2R = 0;
    end
    if P_b2L ~=0
        b2L = OL * b(i)/P_b2L; %The next possible belief when o_t+1 = O_L
    else
        b2L = 0;
    end
    %Floors and ceils of b2R and b2L
    jR_F = floor(b2R * (nb-1)) + 1;
    jR_C = ceil(b2R * (nb-1)) + 1;
    jL_F = floor(b2L * (nb-1)) + 1;
    jL_C = ceil(b2L * (nb-1)) + 1; 
    %Weights on each grid points
    a_R = b2R * (nb - 1) + 1 - jR_F;    
    a_L = b2L * (nb - 1) + 1 - jL_F;
    T(i,jR_F) = T(i,jR_F) + P_b2R * (1 - a_R);
    T(i,jR_C) = T(i,jR_C) + P_b2R * a_R;
   
    T(i,jL_F) = T(i,jL_F) + P_b2L * (1 - a_L);
    T(i,jL_C) = T(i,jL_C) + P_b2L * a_L;
    
    L(i,1) = b(i) * r_p + (1 - b(i)) * r_n;
    L(i,2) = b(i) * r_n + (1 - b(i)) * r_p;
end

F(:,:,1) = eye(nb);
F(:,:,2) = eye(nb);


for k = 1 : (nu-2)/2
    u =  k * 2 + 1;
    F(:,:,u) = F(:,:,u-2) * T;
    F(:,:,u + 1) = F(:,:,u);
    L(:,u) = r_s * k + L(:,1);
    L(:,u+1) = r_s * k + L(:,2);
end



