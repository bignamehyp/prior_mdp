function [ Prop, RT, v, STD] = randomDots(b0, c, len, policy, w_prior)
%Output
%    P: Proportional Positive Choices
%    R: Reaction Time
%    v: decision variable that minicks the LIP firing rate
%Input
%   b0: initial belief
%    c: motion strength
%  len: length of observed sequence
%    d: decision policy


%Interpolate decision boundaries
%policy = interp1(linspace(1,len,length(d)), d,[1:len]);
u =  ( 1+ c )/2;
numTrials = 10000;
P = 0;
N = 0;
iter = 0;
R_P = 0;
R_N = 0;
STD_P = 0;
STD_N = 0;
v = zeros(1,len);
if size(policy, 2) == 1
    policy = policy';
end;
while iter < numTrials
    iter = iter + 1;
    %Generate random observation sequnce o
    o = rand(1,len);
    o(o > u) = 1;
    o(o <=  u) = 0;
    o = 1 - o;
  
    b2 = (cumsum(o)*(1 - w_prior) + b0 * w_prior * (1:len)) ./ (1:len);
    tR = find(b2 > policy);
    tL = find(b2 < 1 - policy);
    
    if isempty(tR)
        tR(1) = len;
    end
    
    if isempty(tL)
        tL(1) = len;
    end
    
    if tR(1) < tL(1)
        P = P + 1;
        R_P = R_P + tR(1);
        STD_P = STD_P + tR(1) * tR(1);
    end    
    
   if tL(1) < tR(1)
        N = N + 1;
        R_N = R_N + tL(1);
        STD_N = STD_N + tL(1) * tL(1);
   end
    
%     nR = 0;
%     nL = 0;
%    for t = 1  : len        
%         nR = nR + o(t);
%         nL = nL + 1 - o(t);
%         d = policy(t+1,:);
%         d(d==0) = [];
%         if d(nR+1) == 2
%             P = P + 1;
%             R_P = R_P + t;
%             STD_P = STD_P + t * t;
%             break;
%         elseif d(nR+1) == 3
%             N = N + 1;
%             R_N = R_N + t;
%             STD_N = STD_N + t * t;
%             break;
%         end
%    end

end
R_P = R_P / P;
R_N = R_N / N;
STD_P = sqrt((STD_P - R_P * R_P * P) / (P - 1));
STD_N = sqrt((STD_N - R_N * R_N * N) / (N - 1));
v = v / iter;
Prop = [P N];
RT = [R_P R_N];
STD = [STD_P STD_N];
end


