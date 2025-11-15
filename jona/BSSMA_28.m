% Backtracking Search Optimization Algorithm (BSA) modified by Jona
% 2023-11-12.
function [best_pos,Convergence_curve] = BSSMA_28(N, Max_FEs, lb, ub, dim, fobj)
disp('BSSMA_28 is tackling your problem.');
tic

%INITIALIZATION
DIM_RATE = 1;  %此变量的取值在原文中并没有给出
Convergence_curve=[];
FEs = 0;
bestFitness = inf;
best_pos = zeros(1,dim);
iter = 1;
AllFitness = inf * ones(N, 1);
if numel(lb)==1, lb=lb*ones(1,dim); ub=ub*ones(1,dim); end
X=initialization(N,dim,ub,lb); % see Eq.1 in [1]
for i = 1 : N
    AllFitness(i) = fobj(X(i, :));
    FEs  = FEs +1;
    if AllFitness(i) < bestFitness
        bestFitness = AllFitness(i);
        best_pos = X(i,:);
    end
end
% historical_X  is swarm-memory of BSA as mentioned in [1].
historical_X = initialization(N, dim, ub, lb);% see Eq.2 in [1]

while FEs < Max_FEs / 3

    %SELECTION-I
    if rand<rand, historical_X=X; end  % see Eq.3 in [1] redefine oldP at the beginning of each iteration
    historical_X=historical_X(randperm(N),:); % see Eq.4 in [1] randomly change the order of the individuals in oldP
    F=get_scale_factor; % see Eq.5 in [1], you can define other F generation strategies 
    map=zeros(N,dim); % see Algorithm-2 in [1]  我看论文中这里是定义的全为1的矩阵呢？       
    if rand<rand
        for i=1:N,  u=randperm(dim); map(i,u(1:ceil(DIM_RATE*rand*dim)))=1; end
    else
        for i=1:N,  map(i,randi(dim))=1; end
    end

    % RECOMBINATION (MUTATION+CROSSOVER)   
    offsprings=X+(map.*F).*(historical_X-X);   % see Eq.5 in [1]  
    %原文中不是X+F*(historical_X-X)吗？对，但这里将变异与交叉合并在一起了，  
    offsprings=BoundaryControl(offsprings,lb,ub); % see Algorithm-3 in [1]

    % SELECTON-II
    for i = 1:N
        fitnessoffsprings(i) = feval(fobj,offsprings(i,:));
        FEs = FEs + 1;
        if fitnessoffsprings(i) < bestFitness
            bestFitness = fitnessoffsprings(i);
            best_pos = offsprings(i,:);
        end
    end
    fitnessoffsprings = reshape(fitnessoffsprings, [1, N]);
    AllFitness = reshape(AllFitness, [1, N]);
    ind=fitnessoffsprings<AllFitness;
    AllFitness(ind)=fitnessoffsprings(ind);
    X(ind,:)=offsprings(ind,:);
    [minimum_fitness, ind]=min(AllFitness);    
    best_pos=X(ind,:);
    bestFitness = minimum_fitness;

    Convergence_curve(iter)=bestFitness;
    iter =iter + 1;
end


% Slime Mold Algorithm modified by Jona 2023-11-10.
% function [best_pos,convergence_curve]=SMA(N,Max_FEs,lb,ub,dim,fobj)

% Initialization
% best_pos=zeros(1,dim);
Destination_fitness=inf;
% AllFitness = inf*ones(N,1);
weight = ones(N,dim);
% X=initialization(N,dim,ub,lb);
% Convergence_curve=[];
% iter=1;  
% lb=ones(1,dim).*lb; 
% ub=ones(1,dim).*ub; 
z=0.03; 

% Main loop
% FEs = 0;
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
                    X(i,j) = best_pos(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    Convergence_curve(iter)=Destination_fitness;
    iter=iter+1;
end
toc
end


function X=BoundaryControl(X,low,up)
[N,dim]=size(X);
for i=1:N
    for j=1:dim                
        k=rand<rand; % you can change boundary-control strategy
        if X(i,j)<low(j)
            if k, X(i,j)=low(j); 
            else X(i,j)=rand*(up(j)-low(j))+low(j); 
            end 
        end        
        if X(i,j)>up(j)
            if k, X(i,j)=up(j);  
            else
                X(i,j)=rand*(up(j)-low(j))+low(j); 
            end 
        end
    end
end
end



function F=get_scale_factor % you can change generation strategy of scale-factor,F    
     F=3*randn; % STANDARD brownian-walk
    % F=4*randg;  % brownian-walk    
    % F=lognrnd(rand,5*rand);  % brownian-walk              
    % F=1/normrnd(0,5);        % pseudo-stable walk (levy-like)
    % F=1./gamrnd(1,0.5);      % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)   
end