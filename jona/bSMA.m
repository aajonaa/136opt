% Slime mold algorithm binary version modified by Jona 20230-10-27.
function [bestFitness,best_pos,convergence_curve, Time]=bSMA(N, Max_iter, dim, A, trn, vald, ~, classifierFhd)
tic
disp('SMA is now tackling your problem')

TFid = 8;

% initialize position
best_pos=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(1, N);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
ub = 1;
lb = 0;
X=initialization(N,dim,ub,lb) > 0.5;
convergence_curve=[];
% it=1;  %Number of iterations
% lb=ones(1,dim).*lb; % lower boundary 
% ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

% Main loop
iter = 1;
while  iter <= Max_iter
    
    %sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
%         Flag4ub=X(i,:)>ub;
%         Flag4lb=X(i,:)<lb;
%         X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = AccSz2(X(i,:), A, trn, vald, classifierFhd);
%         FEs = FEs + 1;
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
    
    a = atanh(-(iter/Max_iter)+1);   %Eq.(2.4)
    b = 1-iter/Max_iter;
    % Update the Position of search agents
    for i=1:N
        if rand<z     %Eq.(2.7)
            for j = 1:dim
                X(i,j) = (ub-lb)*rand+lb;
                X(i, j) = transferFun(X(i, j), X(i, j), TFid);
            end
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(2.2)
            vb = unifrnd(-a,a,1,dim);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,dim);
            for j=1:dim
                r = rand();
                AA = randi([1,N]);  % two positions randomly selected from population
                B = randi([1,N]);
                if r<p    %Eq.(2.1)
                    X(i,j) = best_pos(j)+ vb(j)*(weight(i,j)*X(AA,j)-X(B,j));
                    X(i, j) = transferFun(X(i, j), X(i, j), TFid);
                else
                    X(i,j) = vc(j)*X(i,j);
                    X(i, j) = transferFun(X(i, j), X(i, j), TFid);
                end
            end
        end
    end
    convergence_curve(iter)=Destination_fitness;
    iter=iter+1;
end
Time = toc;
end
