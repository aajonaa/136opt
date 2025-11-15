%% Slime Mold Algorithm modified by Jona 2023-11-20.
% Select the best population from two population.




function [best_pos,Convergence_curve]=BS_SMA_28_2(N,Max_FEs,lb,ub,dim,fobj)
tic
disp('BS_SMA_28_2 is now tackling your problem')

%% initialize position
best_pos=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(1,N);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

oldX = initialization(N, dim, ub, lb);

Convergence_curve=[];
it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter
FEs = 0;

%% Check if solutions go outside the search space and bring them back
for i = 1:N
    Flag4ub=X(i,:)>ub;
    Flag4lb=X(i,:)<lb;
    X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    AllFitness(i) = fobj(X(i, :));
    FEs = FEs + 1;
end

%% Main loop
while  FEs < Max_FEs
    
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);

    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero

    %% calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %% update the best fitness value and best position
    if bestFitness < Destination_fitness
        best_pos=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
    %%  Update the Position of search agents
    a = atanh(-(FEs/Max_FEs)+1);   %Eq.(2.4) a is between [-inf, inf]. 1 - FEs / Max_FEs is between [-1, 1]. 
    b = 1-FEs/Max_FEs;
    for i=1:N
        % Random distribution status.
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        % Exploration status.
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(2.2) p is between [0, 1]. Decreasing from 1 to 0. 1-->0
            vb = unifrnd(-a,a,1,dim);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,dim);
            for j=1:dim
                r = rand();
                A = randi([1,N]);  % two positions randomly selected from population
                B = randi([1,N]);
                % Exploration status.
                if r<p    %Eq.(2.1)
                    X(i,j) = best_pos(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                % Exploitation status. 
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    
    %% Check if updated agents go outside the search space and bring them back
    for i = 1:N
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
    end
    
    %% Selection I
    if rand < rand 
        oldX = X;
    end
    oldX = oldX(randperm(N), :);
    %% Mutation and Crossover
    map = zeros(N, dim);
    if rand < rand
        for i = 1:N
            u = randperm(dim);
            map(i, u(1:ceil(rand * dim))) = 1;
        end
    else
        for i = 1:N
            map(i, randi(dim)) = 1;
        end
    end
    F = 3 * randn();
    offsprings = X + (map .* F) .* (oldX - X);
    offsprings = BoundaryControl(offsprings, lb, ub);
    
    %% Selection II
    fitness_offsprings = inf * ones(1, N);
    for i = 1:N
        fitness_offsprings(i) = fobj(offsprings(i, :));
    end
%     ind = fitness_offsprings < AllFitness;
%     AllFitness(ind) = fitness_offsprings(ind);
%     X(ind, :) = offsprings(ind, :);
    fitness_all = [AllFitness, fitness_offsprings];
    X_all = [X; offsprings];
    [sorted_fitness, ind] = sort(fitness_all);
    AllFitness = sorted_fitness(1:N);
    X = X_all(ind(1:30), :);
    [minimumFitness, index] = min(AllFitness);
    if minimumFitness < Destination_fitness
        Destination_fitness  = minimumFitness;
        best_pos = X(index, :);
    end
    
    Convergence_curve(it)=Destination_fitness;
    it=it+1;
end
toc
end

%% Boundries control
function X=BoundaryControl(X,lb,ub)
[N,dim]=size(X);
for i=1:N
    for j=1:dim                
        k=rand<rand; % you can change boundary-control strategy
        if X(i,j)<lb(j)
            if k, X(i,j)=lb(j); 
            else X(i,j)=rand*(ub(j)-lb(j))+lb(j); 
            end 
        end        
        if X(i,j)>ub(j)
            if k, X(i,j)=ub(j);  
            else
                X(i,j)=rand*(ub(j)-lb(j))+lb(j); 
            end 
        end
    end
end
end
