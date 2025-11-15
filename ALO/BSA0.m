 % Backtracking Search Optimization Algorithm (BSA) modified by Jona
% 2024-3-12.
function [best_pos,Convergence_curve] = BSA0(N, Max_FEs, lb, ub, dim, fobj)
    tic
    
    %INITIALIZATION
    DIM_RATE = 1;  
    Convergence_curve = [];
    FEs = 0;
    bestFitness = inf;
    best_pos = zeros(1,dim);
    iter = 1;
    AllFitness = inf * ones(1, N);
    fitnessoffsprings = inf * ones(1, N);
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
    
    while FEs < Max_FEs
    
        %SELECTION-I
        if rand<rand, historical_X=X; end  % see Eq.3 in [1] redefine oldP at the beginning of each iteration
        historical_X=historical_X(randperm(N),:); % see Eq.4 in [1] randomly change the order of the individuals in oldP
        F=get_scale_factor; % see Eq.5 in [1], you can define other F generation strategies 
        map=zeros(N,dim); % see Algorithm-2 in [1]  
        if rand<rand %0.5
            for i=1:N,  u=randperm(dim); map(i,u(1:ceil(DIM_RATE*rand*dim)))=1; end
        else
            for i=1:N,  map(i,randi(dim))=1; end
        end
    
        % RECOMBINATION (MUTATION+CROSSOVER)   
        offsprings=X+(map.*F).*(historical_X-X);   % see Eq.5 in [1]  
        offsprings=BoundaryControl(offsprings,lb,ub); % see Algorithm-3 in [1]
    
        % SELECTON-II
        for i = 1:N
            fitnessoffsprings(i) = fobj(offsprings(i,:));
            FEs = FEs + 1;
            if fitnessoffsprings(i) < bestFitness
                bestFitness = fitnessoffsprings(i);
                best_pos = offsprings(i,:);
            end
        end
        ind=fitnessoffsprings<AllFitness;
        AllFitness(ind)=fitnessoffsprings(ind);
        X(ind,:)=offsprings(ind,:);
        [minimum_fitness, ind]=min(AllFitness);    
        best_pos=X(ind,:);
        bestFitness = minimum_fitness;
    
        Convergence_curve(iter)=bestFitness;
        iter =iter + 1;
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