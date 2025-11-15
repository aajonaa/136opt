% Backtracking Search Algorithm modified by Jona 2023-10-29.
function  [best_pos,Convergence_curve] = BSA(N,Max_FEs,lb,ub,dim,fobj)
tic

%% INITIALIZATION
DIM_RATE = 1;  
Convergence_curve=[];
FEs = 0;
% best_fitness = 1;
bestFitness = inf;
best_pos = zeros(1,dim);
if numel(lb)==1
    lb=lb*ones(1,dim); 
    ub=ub*ones(1,dim); 
end 
AllFitness = inf * ones(1, N);
X=initialization(N,dim,ub,lb); % see Eq.1 in [1]
% historical_X  is swarm-memory of BSA as mentioned in [1].
historyX=initialization(N,dim,ub,lb);% see Eq.2 in [1]
it = 1;

% Caculate the fitness of the initial population
for i = 1:N
    AllFitness(i) = fobj(X(i, :));
    FEs  = FEs +1;
    if AllFitness(i) < bestFitness
        bestFitness = AllFitness(i);
        best_pos = X(i,:);
    end
end

%% Main loop, include selection mutation and crossover
while FEs < Max_FEs
    %% SELECTION-I
    % Selection mechanism(the operation is part of the algorithm's
    % selection mechanism : based on the random condition)
    if rand < rand
        historyX=X; 
    end  % see Eq.3 in [1] redefine oldP at the beginning of each iteration
    
    % Shuffle the rows of the historical_X matrix, this operation
    % rearranges the order of the solutions in the historical_X matrix.
    % Prevent the algorithm from getting stuck in local optima.
    historyX=historyX(randperm(N),:); % see Eq.4 in [1] randomly change the order of the individuals in oldP
    F=get_scale_factor; % see Eq.5 in [1], you can define other F generation strategies 
    map=zeros(N,dim); % see Algorithm-2 in [1]  
    if rand<rand
        for i=1:N  
            % Permutation the 30 number
            u=randperm(dim); 
            
            % DIM_RATE is a parameter that determines the mutation rate or
            % the propotion of decision variables to be muted.
            % (1:ceil(DIM_RATE*rand*dim) means the vector that should be
            % the index
            % 1:ceil(DIM_RATE*rand*dim) from 1 column to a integer number that randomly selected by the rand 
            map(i,u(1:ceil(DIM_RATE*rand*dim)))=1;
        end
    else
        for i=1:N  
            % randi returns a scalar number that between 1 and dim
            map(i,randi(dim))=1;
        end
    end
    
    %% RECOMBINATION (MUTATION+CROSSOVER)   
    % Get the next generation by the mutation and crossover.
    % Generating a new set of candidate solutions by applying mutation and
    % recombination operation to the current popolation(X) and the
    % historical population(historical_X)
    % Map is a binary matrix that determines which components of the
    % solution should undergo mutation.
    % F is a scaling factor that controls the magnitude of the mutation and
    % direction of the mutation.
    newX=X+(map.*F).*(historyX-X);   % see Eq.5 in [1]  
    newX=BoundaryControl(newX,lb,ub); % see Algorithm-3 in [1]
    
    %% SELECTON-II
    % Select the best individual based on the next generation's fitness to
    % the problem.
    for i = 1:N
        newAllFitness(i) = fobj(newX(i, :));
        FEs = FEs + 1;
    end
    
    % Create a logical index or mask based on a comparison between two
    % array. To determine whether the fitness values of newly generated
    % candidate solutions are better than the fitness values of the current
    % population.
    ind = newAllFitness < AllFitness;
    AllFitness(ind) = newAllFitness(ind);
    % We select better solution of offspring to assign it to the origin
    % population.
    X(ind, :)=newX(ind, :);
    % Get the global minimum(the minimal fitness) and the index of the best
    % solution.
    [globalminimum, ind] = min(AllFitness);
    % Get the best solution(the best individual: the solution row).
    globalminimizer = X(ind, :);
    bestFitness = globalminimum;
    best_pos = globalminimizer;
    Convergence_curve(it) = bestFitness;
    it =it + 1;
end
Time = toc;
disp(['The BSA run for ', num2str(Time), 'Second']);
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

%% F factor
function F=get_scale_factor % you can change generation strategy of scale-factor,F    
     F=3*randn; % STANDARD brownian-walk
    % F=4*randg;  % brownian-walk    
    % F=lognrnd(rand,5*rand);  % brownian-walk              
    % F=1/normrnd(0,5);        % pseudo-stable walk (levy-like)
    % F=1./gamrnd(1,0.5);      % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)   
end