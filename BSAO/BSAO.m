%% The original AO.

function [Leader_pos,Convergence_curve]=AO(N,Max_FEs,lb,ub,dim,fobj)
    % function [bestfitness,Convergence_curve]=AO0(N,MaxFEs,lb,ub,dim,fobj) % For myself
    %初始化参数
    FEs=0;
    it=1;
    Fitnorm=zeros(1,N);
    Convergence_curve=[];
    %% 种群的初始化
    % 初始化一个个体
    X=initialization(N,dim,ub,lb);
    %计算初始种群的适应度值
    for i=1:N
        Fitness(i)=fobj(X(i,:));
        FEs=FEs+1;
    end
    % 为初始种群排序成新种群，找出最优个体，并记录
    [fmin,x]=min(Fitness);
    %
    New_pop=zeros(N,dim);
    best=X(x,:);
    bestfitness=fmin;
    %% 主循环
    while FEs<=2/3 * Max_FEs
        
        K= 1-((FEs)^(1/6)/(Max_FEs)^(1/6));
        %E 是一个随迭代递减的数值，它可以成为控制后期局部开发的权重，调节全局搜索和局部开发
        E =1*exp(-4*(FEs/Max_FEs));
    
        for i=1: N
            
            Fitnorm(i)= (Fitness(i)-min(Fitness))/(max(Fitness)-min(Fitness));
            for j=1:dim
                %攻击阶段  
                if rand<K %(Comprehensive elimination phase
                    if rand<0.5
                        New_pop(i,j) = X(i,j)+E.*X(i,j)*(-1)^FEs; %% 1
                    else
                        New_pop(i,j) = X(i,j)+E.*best(j)*(-1)^FEs; %% 2
                    end
                else
                    New_pop(i,j)=X(i,j);
                end
                
                if rand<Fitnorm(i) %(Local clearance phase
                    A=randperm(N);
                    beta=(rand/2)+0.1;
                    New_pop(i,j)=X(A(3),j)+beta.*(X(A(1),j)-X(A(2),j)); %% 3
                end
            end
            
            New_pop(i,:)=Mutation(New_pop(i,:),X(i,:),best,dim); %% 4 %(Post-consolidation phase
            %% 边界收束 ：重新给药，并不是所有的二氧青蒿素都会寻找到疟原虫，一些会因为自身代谢排出人体，需要再次给药、
            New_pop(i,:)=Transborder_reset(New_pop(i,:),ub,lb,dim,best); 
            % 计算适应度值
            tFitness=fobj(New_pop(i,:));
            FEs=FEs+1;
            %初步更新适应度值
            if tFitness<Fitness(i)
                X(i,:)= New_pop(i,:);
                Fitness(i)=tFitness;
            end
        end
        [fmin,x]=min(Fitness);
        if fmin<bestfitness
            best=X(x,:);
            bestfitness=fmin;
        end
        % 保存每次迭代的最佳函数值
        Convergence_curve(it)=bestfitness;
        Leader_pos=best;
        bestfitness = min(Fitness);
        it=it+1;
    end
% end

% Backtracking Search Algorithm modified by Jona 2023-10-29.
% function  [best_pos,convergence_curve] = BSA(N,Max_FEs,lb,ub,dim,fobj)
tic

%% INITIALIZATION
DIM_RATE = 1;  
% convergence_curve=[];
% FEs = 0;
% best_fitness = 1;
% bestfitness = inf;
% Leader_pos = zeros(1,dim);
if numel(lb)==1
    lb=lb*ones(1,dim); 
    ub=ub*ones(1,dim); 
end 
fitnessX = Fitness;
% X=initialization(N,dim,ub,lb); % see Eq.1 in [1]
% historical_X  is swarm-memory of BSA as mentioned in [1].
% historical_X=initialization(N,dim,ub,lb);% see Eq.2 in [1]
historical_X=New_pop;% see Eq.2 in [1]
% it = 1;

% Caculate the fitness of the initial population
% for i = 1:N
%     fitnessX(i) = fobj(X(i, :));
%     FEs  = FEs +1;
%     if fitnessX(i) < bestfitness
%         bestfitness = fitnessX(i);
%         Leader_pos = X(i,:);
%     end
% end

%% Main loop, include selection mutation and crossover
while FEs < Max_FEs
    %% SELECTION-I
    % Selection mechanism(the operation is part of the algorithm's
    % selection mechanism : based on the random condition)
    if rand < rand
        historical_X=X; 
    end  % see Eq.3 in [1] redefine oldP at the beginning of each iteration
    
    % Shuffle the rows of the historical_X matrix, this operation
    % rearranges the order of the solutions in the historical_X matrix.
    % Prevent the algorithm from getting stuck in local optima.
    historical_X=historical_X(randperm(N),:); % see Eq.4 in [1] randomly change the order of the individuals in oldP
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
    offsprings=X+(map.*F).*(historical_X-X);   % see Eq.5 in [1]  
    offsprings=BoundaryControl(offsprings,lb,ub); % see Algorithm-3 in [1]
    
    %% SELECTON-II
    % Select the best individual based on the next generation's fitness to
    % the problem.
    for i = 1:N
        fitness_offsprings(i) = fobj(offsprings(i, :));
        FEs = FEs + 1;
    end
    
    % Create a logical index or mask based on a comparison between two
    % array. To determine whether the fitness values of newly generated
    % candidate solutions are better than the fitness values of the current
    % population.
    ind = fitness_offsprings < fitnessX;
    fitnessX(ind) = fitness_offsprings(ind);
    % We select better solution of offspring to assign it to the origin
    % population.
    X(ind, :)=offsprings(ind, :);
    % Get the global minimum(the minimal fitness) and the index of the best
    % solution.
    [globalminimum, ind] = min(fitnessX);
    % Get the best solution(the best individual: the solution row).
    globalminimizer = X(ind, :);
    bestfitness = globalminimum;
    Leader_pos = globalminimizer;
    Convergence_curve(it) = bestfitness;
    it =it + 1;
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

%% F factor
function F=get_scale_factor % you can change generation strategy of scale-factor,F    
     F=3*randn; % STANDARD brownian-walk
    % F=4*randg;  % brownian-walk    
    % F=lognrnd(rand,5*rand);  % brownian-walk              
    % F=1/normrnd(0,5);        % pseudo-stable walk (levy-like)
    % F=1./gamrnd(1,0.5);      % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)   
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

function z=Transborder_reset(z,ub,lb,dim,best)
    for j=1:dim
        if z(j)>ub || z(j)<lb
            
            z(j)=best(j);
            
        end
    end
end