%___________________________________________________________________%
%  Ant Lion Optimizer (ALO) source codes demo version 1.0           %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%   S. Mirjalili, The Ant Lion Optimizer                            %
%   Advances in Engineering Software , in press,2015                %
%   DOI: http://dx.doi.org/10.1016/j.advengsoft.2015.01.010         %
%                                                                   %
%___________________________________________________________________%

% This function creates random walks

function [RWs]=Random_walk_around_antlion(Dim,Max_FEs,lb, ub,antlion,FEs)
if size(lb,1) ==1 && size(lb,2)==1 %Check if the bounds are scalar
    lb=ones(1,Dim)*lb;
    ub=ones(1,Dim)*ub;
end

if size(lb,1) > size(lb,2) %Check if boundary vectors are horizontal or vertical
    lb=lb';
    ub=ub';
end

I=1; % I is the ratio in Equations (2.10) and (2.11)

if FEs>Max_FEs/10
    I=1+100*(FEs/Max_FEs);
end

if FEs>Max_FEs/2
    I=1+1000*(FEs/Max_FEs);
end

if FEs>Max_FEs*(3/4)
    I=1+10000*(FEs/Max_FEs);
end

if FEs>Max_FEs*(0.9)
    I=1+100000*(FEs/Max_FEs);
end

if FEs>Max_FEs*(0.95)
    I=1+1000000*(FEs/Max_FEs);
end


% Dicrease boundaries to converge towards antlion
lb=lb/(I); % Equation (2.10) in the paper
ub=ub/(I); % Equation (2.11) in the paper

% Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
if rand<0.5
    lb=lb+antlion; % Equation (2.8) in the paper
else
    lb=-lb+antlion;
end

if rand>=0.5
    ub=ub+antlion; % Equation (2.9) in the paper
else
    ub=-ub+antlion;
end

% This function creates n random walks and normalize accroding to lb and ub
% vectors 
RWs_iter = zeros(1, Dim);
for i=1:Dim
    X = [0 cumsum(2* (rand(Max_FEs,1)>0.5) -1)']; % Equation (2.1) in the paper
    %[a b]--->[c d]
    a=min(X);
    b=max(X);
    c=lb(i);
    d=ub(i);      
%     X_norm=((X-a).*(d-c))./(b-a)+c; % Equation (2.7) in the paper
    X_norm = ((X(end)- min(X)) *(ub(i) - lb(i)) / (max(X) - min(X))) + lb(i);
%     RWs(:,i)=X_norm;
    RWs_iter(i) = X_norm;
end
RWs = RWs_iter;

