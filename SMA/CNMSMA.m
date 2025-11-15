function [bestPositions,Convergence_curve]=CNMSMA(N,Max_FEs,lb,ub,dim,fobj)
disp('SMA_FEs is now tackling your problem')
t = cputime;
% initialize position
bestPositions=zeros(1,dim);
Destination_fitness=1;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
Convergence_curve=[];
FEs=0;   %Number of current evaluations
MaxFEs=Max_FEs;  %Maximum number of evaluations
fitcount = 0;
it=1;   %Number of iterations
lb=ones(1,dim).*lb; % lower boundary
ub=ones(1,dim).*ub; % upper boundary
z=0.045; % parameter0.03

% Main loop
while  FEs < MaxFEs
    'CNMSMA'
    %sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
        FEs = FEs+1;
        if AllFitness(i) < Destination_fitness 
            Destination_fitness = AllFitness(i);
            bestPositions = X(i,:);
            INV=1;
        end
        fitcount = fitcount + 1;
        Convergence_curve(1,fitcount) = Destination_fitness;
    end
    
        %%%混沌局部搜索***%%%、
        setCan = (Max_FEs-FEs+1)/Max_FEs;
        x = rand();
        if(x~=0.25&&x~=0.5&&x~=0.75)
            ch(1) = x;
        end
        for ant=1:(N)
            ch(ant+1)=4*ch(ant)*(1-ch(ant));
            CH(ant,:) = lb+ch(ant)*(ub-lb);    %ub大
            V = (1-setCan)*bestPositions+setCan*CH(ant);
            % ObjValV=feval(objfun,V);         %计算函数值
            [FitnessV]=fobj(V);%计算适应度值
            if FitnessV<fobj(X(i,:))
                X(i,:)=V;
            end
            if (FitnessV<Destination_fitness)
                Destination_fitness = FitnessV;
                bestPositions = V;
                break;
            end
        end
   

    
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);

    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero

    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)    %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
        INV=1;
    end
    
    %% NMs
    options = optimset('MaxFunEvals', floor(MaxFEs * 0.1));
    [x, fval, ~, output]  = fminsearchbnd(fobj,bestPositions,lb,ub,options);
    if fval < Destination_fitness
        Destination_fitness = fval;
        bestPositions = x;
        INV=1;
    end
    FEs = FEs + output.funcCount;
    Convergence_curve(:,fitcount + 1:fitcount + output.funcCount) = Destination_fitness;
    fitcount = fitcount + output.funcCount;
    %%    
    
    a = atanh(-(FEs/Max_FEs)+1); %Eq.(2.4)
    b = 1-FEs/Max_FEs; 
    % Update the Position of search agents
    for i=1:N
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness)); %Eq.(2.2)
            vb = unifrnd(-a,a,1,dim);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,dim);
            for j=1:dim
                r = rand();
                A = randi([1,N]);    % two positions randomly selected from population
                B = randi([1,N]);
                if r<p          %Eq.(2.1)
                    X(i,j) = bestPositions(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    
    
end

Convergence_curve=Convergence_curve(:,(fitcount-MaxFEs+1):end);
time = cputime - t;
end
