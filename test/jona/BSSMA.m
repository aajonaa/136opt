% Backtracking Search Algorithm combined Slime Mould Algorithm reviewed
% second time at third week by Jona 2023-11-2.

function  [best_pos,Convergence_curve] = BSSMA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
tic
disp('BSSMA is tackling your problem');
DIM_RATE = 1;  %此变量的取值在原文中并没有给出
Convergence_curve=[];
FEs = 0;
% best_fitness = 1;
best_fitness = inf;
best_pos = zeros(1,dim);
it = 1;
%INITIALIZATION
if numel(lb)==1, lb=lb*ones(1,dim); ub=ub*ones(1,dim); end % this line must be adapted to your problem
%numel:Number of elements in an array or subscripted array expression.
X=initialization(SearchAgents_no,dim,ub,lb); % see Eq.1 in [1]
% historical_X  is swarm-memory of BSA as mentioned in [1].
historical_X=initialization(SearchAgents_no,dim,ub,lb);% see Eq.2 in [1]



% fitnessX=feval(fobj,X);
% Caculate the fitness of the initial population
for i = 1 : SearchAgents_no
    fitnessX(i) = feval(fobj,X(i,:));
    FEs  = FEs +1;
    if fitnessX(i) < best_fitness
        best_fitness = fitnessX(i);
        best_pos = X(i,:);
    end
end



% Main loop, include selection mutation and crossover
% ------------------------------------------------------------------------------------------ 
while FEs < Max_iter / 3
    %SELECTION-I
    % Selection mechanism(the operation is part of the algorithm's
    % selection mechanism : based on the random condition)
    if rand<rand, historical_X=X; end  % see Eq.3 in [1] redefine oldP at the beginning of each iteration
    
    % Shuffle the rows of the historical_X matrix, this operation
    % rearranges the order of the solutions in the historical_X matrix.
    % Prevent the algorithm from getting stuck in local optima.
    historical_X=historical_X(randperm(SearchAgents_no),:); % see Eq.4 in [1] randomly change the order of the individuals in oldP
    F=get_scale_factor; % see Eq.5 in [1], you can define other F generation strategies 
    map=zeros(SearchAgents_no,dim); % see Algorithm-2 in [1]  我看论文中这里是定义的全为1的矩阵呢？       
    if rand<rand
        for i=1:SearchAgents_no  
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
        for i=1:SearchAgents_no  
            % randi returns a scalar number that between 1 and dim
            map(i,randi(dim))=1;
        end
    end
    
    % RECOMBINATION (MUTATION+CROSSOVER)   
    % Get the next generation by the mutation and crossover.
    % Generating a new set of candidate solutions by applying mutation and
    % recombination operation to the current popolation(X) and the
    % historical population(historical_X)
    % Map is a binary matrix that determines which components of the
    % solution should undergo mutation.
    % F is a scaling factor that controls the magnitude of the mutation and
    % direction of the mutation.
    offsprings=X+(map.*F).*(historical_X-X);   % see Eq.5 in [1]  
    %原文中不是X+F*(historical_X-X)吗？对，但这里将变异与交叉合并在一起了，  
    
    offsprings=BoundaryControl(offsprings,lb,ub); % see Algorithm-3 in [1]
    
    % SELECTON-II
    % Select the best individual based on the next generation's fitness to
    % the problem.
    for i = 1:SearchAgents_no
        fitnessoffsprings(i) = feval(fobj,offsprings(i,:));
        FEs = FEs + 1;
        if fitnessoffsprings(i) < best_fitness
            best_fitness = fitnessoffsprings(i);
            best_pos = offsprings(i,:);
        end
    end
    
    % Create a logical index or mask based on a comparison between two
    % array. To determine whether the fitness values of newly generated
    % candidate solutions are better than the fitness values of the current
    % population.
    ind=fitnessoffsprings<fitnessX;
    fitnessX(ind)=fitnessoffsprings(ind); %将fitnessoffspring中适应度值小的赋值给fitnessX中相应位置
    % We select better solution of offspring to assign it to the origin
    % population.
    X(ind,:)=offsprings(ind,:);
    % Get the global minimum(the minimal fitness) and the index of the best
    % solution.
    [globalminimum,ind]=min(fitnessX);
    % Get the best solution(the best individual: the solution row).
    globalminimizer=X(ind,:);
    best_fitness = globalminimum;
    best_pos =globalminimizer;
    Convergence_curve(it)=best_fitness;
    it =it + 1;
end

% Slime Mold Algorithem
% function [best_pos,Convergence_curve]=SMA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize position
% best_pos=zeros(1,dim);
best_pos =globalminimizer;
% best_fitness=inf;%change this to -inf for maximization problems
best_fitness = globalminimum;
% AllFitness = inf*ones(SearchAgents_no,1);%record the fitness of all slime mold
AllFitness = fitnessX;
weight = ones(SearchAgents_no,dim);%fitness weight of each slime mold
weight1 = ones(SearchAgents_no,dim);
%Initialize the set of random solutions
% X=initialization(SearchAgents_no,dim,ub,lb);

% X = X;

% Convergence_curve=[];

% Convergence_curve=Convergence_curve;

% it=1;
% FEs = Max_iter / 3 + 1;

% Main loop
while  FEs < Max_iter;
    
    
    %sort the fitness
    
    for i=1:size(X,1)
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        FEs = FEs + 1;
        AllFitness(i) = fobj(X(i,:));
    end
    [SmellOrder,SmellIndex] = sort(AllFitness);
    worstFitness = SmellOrder(size(X,1));
    bestFitness = SmellOrder(1);
    

    %calculate the fitness weight of each slime mold
    for i=1:size(X,1)
        for j=2:size(X,2)
            if i<=(size(X,1)/2)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(bestFitness-worstFitness)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(bestFitness-worstFitness)+1);
            end
            if SmellOrder(i)>mean(AllFitness)
                weight1(i,:) = 1+(worstFitness-SmellOrder(i))/(worstFitness-bestFitness);
            else
                weight1(i,:) = 1.000001-(worstFitness-SmellOrder(i))/(worstFitness-bestFitness);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < best_fitness
        best_pos=X(SmellIndex(1),:);
        best_fitness = bestFitness;
    end
    
    a=2*(1-FEs/Max_iter); % a decreases linearly fron 2 to 0
    % Update the Position of search agents
    for i=1:size(X,1)
        if rand<0.03
            X(i,:) = (ub-lb)*rand+lb;
        else
            for j=1:size(X,2)
                r = rand();
                nd = 2*a*r-a;%[-a,a]
                D = nd*abs(weight(i,j)*X(i,j)-best_pos(j));
                C=2*rand();  %[0,2]
                if r>0.1
                    X(i,j) = best_pos(j)-D;
                else
                    X(i,j) = weight1(i,j)*(best_pos(j)-nd*abs(C*best_pos(j)-X(i,j))); 
                end
            end
        end
    end
    
    Convergence_curve(it)=best_fitness;
    it=it+1;
end

end

function X=BoundaryControl(X,lb,ub)
[SearchAgents_no,dim]=size(X);
for i=1:SearchAgents_no
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
toc
end



function F=get_scale_factor % you can change generation strategy of scale-factor,F    
     F=3*randn; % STANDARD brownian-walk
    % F=4*randg;  % brownian-walk    
    % F=lognrnd(rand,5*rand);  % brownian-walk              
    % F=1/normrnd(0,5);        % pseudo-stable walk (levy-like)
    % F=1./gamrnd(1,0.5);      % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)   
% return
end
