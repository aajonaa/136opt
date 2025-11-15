%% Multi-strategy AO developed by Jona 2024-6.16.

function [best_pos,Convergence_curve]=MSAO(N,Max_FEs,lb,ub,dim,fobj) 
% function [bestFitness,Convergence_curve]=AO(N,Max_FEs,lb,ub,dim,fobj) % For myself

FEs=0;
it=1;
FitNorm= inf * ones(N, 1);
newFitAll = inf * ones(N, 1);
FitTrack = inf * ones(N, 3);
FitAll = inf * ones(N, 1);
Convergence_curve=[];

X=initialization(N,dim,ub,lb);

for i=1:N
    FitAll(i)=fobj(X(i,:));
    FEs=FEs+1;
end

% Init the fitTrack
for i = 1:3
    FitTrack(:, i) = FitAll;
end

[fit_min,minfit_idx]=min(FitAll);
%
newX=zeros(N,dim);
best_position=X(minfit_idx,:);
bestFitness=fit_min;

while FEs<=Max_FEs
    
    K= 1-((FEs)^(1/3)/(Max_FEs)^(1/3));
    E =1*exp(-4*(FEs/Max_FEs));
    
    %% Traveling distance rate, the original constant is 1/6
    TDR = K;

    for i=1: N
        
        FitNorm(i)= (FitAll(i)-min(FitAll))/(max(FitAll)-min(FitAll));
        for j=1:dim
 
            if rand<K
%                 if rand<rand
%                     newX(i,j) = X(i,j)+E * X(i,j)*(-1)^FEs; %% 1
%                 else
                    %% Third part change: Group guide 
%                     newX(i,j) = X(i,j)+E * best_position(j)*(-1)^FEs; %% 2
                    [~, index] = sort(FitAll);
                    NumOfAlpha = floor(N/3);
                    tempIdx = zeros(1, NumOfAlpha);
                    for Idx = 1:NumOfAlpha
                        tempIdx(Idx) = index(Idx);
                    end
                    choose = randperm(NumOfAlpha);
                    
                    newX(i, j) = X(i, j) - E * (X(randi(dim), j) - X(i, j)) ...
                        + E * (X(tempIdx(choose(1)), j) - X(i, j)); %% 2
%                     newX(i, j) = X(i, j) + E * (best_position(j) - X(randi(dim), j)) * (-1)^FEs; %% 2
%                 end
            else
                %% First part change: Add the sine-learning mechanism (SLM) based on von Mises distribution (VMD)
                theta = pi/2 * min(FitNorm)/FitNorm(i);
%                 newX(i,j)=X(i,j); %% 3
                newX(i,j) = X(i,j) + sin(theta) * best_position(j); %% 3
            end
        end       
    end
    
    % Boundary control
    newX(i,:)=Transborder_reset(newX(i,:), lb, ub, dim, best_position); %% Boundary control
    % Update the FitAll based on greedy selection
    for i = 1:N
        tempFit = fobj(newX(i, :));
        FEs = FEs + 1;
        if tempFit<FitAll(i)
            X(i,:)= newX(i,:);  %% Greedy selection
            FitAll(i)=tempFit;
        end
    end
    if min(FitAll) < bestFitness
        [~, idx] = min(FitAll);
        best_position = X(idx, :);
    end
    
    for i = 1:N
        %% Second part change: Two stage mutation
        FitTrackNorm = (FitTrack - min(FitTrack)) ./ (max(FitTrack) - min(FitTrack));
        alpha = (FitTrackNorm(i, 2) - FitTrackNorm(i, 1)) / ...
            (sum(FitTrackNorm(:, 2) - FitTrackNorm(:, 1)) + eps)...
            + 0.2;
        factor = 1/2;
        if FEs < factor * Max_FEs
            newX(i,:)=Mutation(newX(i,:),X(i,:),best_position,dim); %% 5
        else
            newX(i, :) = best_position + alpha * TDR * randsrc(1, dim)...
                .* ((ub -lb) .* rand(1, dim) + lb); %% 5.2
        end
        newX(i,:)=Transborder_reset(newX(i,:), lb, ub, dim, best_position); %% Boundary control
        newFitAll(i)=fobj(newX(i,:));
        FEs=FEs+1;
    end
    
    % Track the generation's fitness
    FitTrack(:, 3) = FitTrack(:, 2);
    FitTrack(:, 2) = FitTrack(:, 1);
    FitTrack(:, 1) = newFitAll;
    
    % Update the FitAll based on greedy selection
    for i = 1:N
        if newFitAll(i)<FitAll(i)
            X(i,:)= newX(i,:);  %% Greedy selection
            FitAll(i)=newFitAll(i);
        end
    end
    
    [fit_min,minfit_idx]=min(FitAll);
    if fit_min<bestFitness
        best_position=X(minfit_idx,:);
        bestFitness=fit_min;
    end
    Convergence_curve(it)=bestFitness;
    best_pos=best_position;
    bestFitness = min(FitAll);
    it=it+1;
end
end

function z=Mutation(z,x,b,dim)
    for j=1:dim
        if rand<0.05
            z(j)=x(j);
        end
        if rand<0.2
            z(j)=b(j);
        end
    end
end

function z=Transborder_reset(z,lb, ub, dim, best)
    for j=1:dim
        if z(j)>ub || z(j)<lb || isnan(z(j))
            z(j)=best(j);
        end
    end
end