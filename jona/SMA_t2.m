% Slime Mold Algorithm modified by Jona 2023-11-16.
% 1. sine map to get better initilized population.
function [best_pos,convergence_curve]=SMA_t2(N,Max_FEs,lb,ub,dim,fobj)
tic
disp('SMA_t2 is now tackling your problem')

% initialize position
best_pos=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
% X=initialization(N,dim,ub,lb);

% %% sine map
% for i = 1:N
%     if i == 1
%         for j = 1:dim
%             X(i, j) = (ub - lb) * rand + lb;
%         end
%     else
%         a = rand();
%         for j = 1:dim
%             X(i, j) = (4 * sin(pi * X(i - 1, j))) / (4 * a);
%         end
%     end
% end

% %% logistic tent map
% % mu = 4;
% % gamma = 0.5;
% for i = 1:N
%     mu = 4 * rand();
%     gamma = 0.5 * rand();
%     for j = 1:dim
%         x = lb + (ub - lb) * rand();
% %         X(i, j) = mu * x * (1 - x) * gamma * min(x, 1-x);
%         X(i, j) = mu * (x - lb) / (ub - lb) * (1 - (x - lb) / (ub - lb)) * gamma * min((x - lb) / (ub - lb), 1 - (x - lb) / (ub - lb));
%     end
% end

% %% logistic map
% % mu = 4;
% for i = 1:N
%     if i == 1
%         for j = 1:dim
%             X(i, j) = (ub - lb) * rand() + lb;
%         end
%     else
%         mu = 4 * rand();
%         for j = 1:dim
%             X(i, j) = mu * (X(i-1, j) - lb) / (ub - lb) * (1 - (X(i-1, j) - lb) / (ub - lb));
%         end
%     end
% end

%% tent map
% gamma = 0.5;
for i = 1:N
    if i == 1
        for j = 1:dim
            X(i, j) = (ub - lb) * rand() +lb;
        end
    else
        gamma = 4 * rand();
        for j = 1:dim
            X(i, j) = 4 / gamma * min((X(i - 1, j) - lb) / (ub - lb), 1 - (X(i-1, j) - lb) / (ub - lb));
        end
    end
end

% %% logistic tent map
% % mu = 4;
% % gamma = 0.5;
% for i = 1:N
%     if i == 1
%         for j = 1:dim
%             X(i, j) = (ub - lb) * rand() + lb;
%         end
%     else
%         mu = 4 * rand();
%         gamma = 0.5 * rand();
%         alpha = 4 * rand();
%         for j = 1:dim
% %         X(i, j) = mu * x * (1 - x) * gamma * min(x, 1-x);
%             X(i, j) = mu * (X(i-1, j) - lb) / (ub - lb) * (1 - (X(i - 1, j) - lb) / (ub - lb)) * gamma * min((X(i-1, j) - lb) / (ub - lb), 1 - (X(i-1, j) - lb) / (ub - lb)) * 4 / alpha;
%         end
%     end
% end
% %%

convergence_curve=[];
it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

% Main loop
FEs = 0;
while  FEs <= Max_FEs
    
    %sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
        FEs = FEs + 1;
    end
    
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);

    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero

    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        best_pos=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
    a = atanh(-(FEs/Max_FEs)+1);   %Eq.(2.4)
    b = 1-FEs/Max_FEs;
    % Update the Position of search agents
    for i=1:N
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(2.2)
            vb = unifrnd(-a,a,1,dim);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,dim);
            for j=1:dim
                r = rand();
                A = randi([1,N]);  % two positions randomly selected from population
                B = randi([1,N]);
                if r<p    %Eq.(2.1)
                    X(i,j) = best_pos(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    convergence_curve(it)=Destination_fitness;
    it=it+1;
end
toc
end
