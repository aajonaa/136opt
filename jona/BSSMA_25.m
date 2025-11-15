% Slime Mold Algorithm modified by Jona 2023-10-27.
% Modified by Jona at 2023-11-9, I use the mutation and crossover at the
% end of the search agent update phase to find better solution though it
% may cause the premature result and preconvergence and fall into local
% optimal. But I should have a try at least and then I will add another
% mechanism to improve the exploration phase.
function [best_pos,convergence_curve]=BSSMA_25(N,Max_FEs,lb,ub,dim,fobj)
tic
disp('BSSMA_25 is now tackling your problem')

% initialize position
best_pos=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
convergence_curve=[];
it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

% Main loop
FEs = 0;
while  FEs < Max_FEs
    
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
                    % Exploration phase, extensive search.
                    X(i,j) = best_pos(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    % Exploitation phase, local search.
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    
    bs_X = X(randperm(N), :);
    F = 3 * randn;
    map = zeros(N, dim);
    for i = 1:30
        u = randperm(dim);
        map(i, u(1:ceil(rand * dim))) = 1;
    end
    offsprings = X + (map .* F) .* (bs_X - X);
    offsprings = BoundaryControl(offsprings, lb, ub);
    fitness_offsprings = inf * ones(1, N);
    for i = 1:N
        fitness_offsprings(i) = fobj(offsprings(i, :));
        FEs = FEs + 1;
    end
    fitness_offsprings = reshape(fitness_offsprings, [1, N]);
    AllFitness = reshape(AllFitness, [1, N]);
    ind = fitness_offsprings < AllFitness;
    AllFitness(ind) = fitness_offsprings(ind);
    X(ind, :) = offsprings(ind, :);
    
    
    convergence_curve(it)=Destination_fitness;
    it=it+1;
end
toc
end


function X=BoundaryControl(X,lb,ub)
[N,dim]=size(X);
for i=1:N
    for j=1:dim                
        k=rand<rand; % you can change boundary-control strategy
        if X(i,j)<lb(j)
            if k, X(i,j)=lb(j); 
            else
                X(i,j)=rand*(ub(j)-lb(j))+lb(j);
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
