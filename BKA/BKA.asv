%________________________________________________________ ________________%
%  Black-winged Kite Algorithm (BKA) source codes                         %
%                                                                         %
%  Developed in MATLAB R2022b                                             %
%                                                                         %
%  Author and programmer:                                                 %
%  Black-winged Kite Algorithm: A nature-inspired meta-heuristic for
%              Solving benchmark functions and Engineering problems                                                                       %
%  e-Mail:                                                                %
%  Artificial Intelligence Review                                                                      %                                                                        %
%  DOI:                                                                   %
%                                                                         %
%_________________________________________________________________________%
%%
%%  Black-winged Kite Algorithm
function [best_pos,Convergence_curve]=BKA(N,MaxFEs,lb,ub,dim,fobj)
%% ----------------Initialize the locations of Blue Sheep------------------%
 X=initialization(N,dim,ub,lb);% Initial population

FEs = 0;
it = 1;

for i =1:N
    AllFitness(i)=fobj(X(i,:));
    FEs = FEs + 1;
end
Convergence_curve=[];
%% -------------------Start iteration------------------------------------%
% for t=1:MaxFEs
while FEs <= MaxFEs
    [~,sorted_indexes]=sort(AllFitness);
    Best_pos=X(sorted_indexes(1),:);
    bestFitness = AllFitness(sorted_indexes(1));
   
%% -------------------Attacking behavior-------------------%
    for i=1:N
        
        n=0.05*exp(-2*(FEs/MaxFEs)^2);
        if p<r
            newX(i,:)=X(i,:)+n.*(1+sin(r))*X(i,:);%%1
        else
            newX(i,:)= X(i,:).*(n*(2*rand(1,dim)-1)+1);%%2
        end
        newX(i,:) = max(newX(i,:),lb);newX(i,:) = min(newX(i,:),ub);%%Boundary checking
%% ------------ Select the optimal fitness value--------------%
        
        newAllFitness(i)=fobj(newX(i,:));
        FEs = FEs + 1;
        if(newAllFitness(i)<AllFitness(i))
            X(i,:) = newX(i,:);
            AllFitness(i) = newAllFitness(i);
        end
%% -------------------Migration behavior-------------------%
 
        m=2*sin(r+pi/2);
        s = randi([1,30],1);
        tempFitness=AllFitness(s);
        ori_value = rand(1,dim);cauchy_value = tan((ori_value-0.5)*pi);
        if AllFitness(i)< tempFitness
            newX(i,:)=X(i,:)+cauchy_value(:,dim).* (X(i,:)-Best_pos);%%3
        else
            newX(i,:)=X(i,:)+cauchy_value(:,dim).* (Best_pos-m.*X(i,:));%%4
        end
        newX(i,:) = max(newX(i,:),lb);newX(i,:) = min(newX(i,:),ub); %%Boundary checking
%% --------------  Select the optimal fitness value---------%
        newAllFitness(i)=fobj(newX(i,:));
        FEs = FEs + 1;
        if(newAllFitness(i)<AllFitness(i))
            X(i,:) = newX(i,:);
            AllFitness(i) = newAllFitness(i);
        end
    end
    %% -------Update the optimal Black-winged Kite----------%
    if(AllFitness<bestFitness)
        bestFitness=AllFitness(i);
        best_pos=X(i,:);
    else
        bestFitness=bestFitness;
        best_pos=Best_pos;
    end
    Convergence_curve(it)=bestFitness;
    it = it + 1;
end
end
