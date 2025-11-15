function [bestPositions,Convergence_curve]=HSMA_WOA(N,Max_FEs,lb,ub,dim,fobj)
disp('HSMA_WOA is now tackling your problem')
% t = cputime;
% initialize position
bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
Convergence_curve=[];
FEs=0;   %Number of current evaluations
MaxFEs=Max_FEs;  %Maximum number of evaluations
% fitcount = 0;
CI=12000;   %Number of iterations
lb=ones(1,dim).*lb; % lower boundary
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter
it=1;
% Main loop
while  FEs < MaxFEs
    'HSMA_WOA'
    if FEs<CI
        
        X=WOA_S(N,X,bestPositions,FEs,MaxFEs);
        
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
            end
        end
        
        continue;
        
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
    end
    
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
        end
    end    
    Convergence_curve(it) = Destination_fitness;
    it = it + 1;
end

% Convergence_curve=Convergence_curve(:,(fitcount-MaxFEs+1):end);
% time = cputime - t;
end

function X=WOA_S(N,X,bestPositions,FEs,MaxFEs)
    
    a=2-FEs*((2)/MaxFEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+FEs*((-1)/MaxFEs);
    
    
    % Update the Position of search agents 
    for i=1:size(X,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(X,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(N*rand()+1);
                    X_rand = X(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-X(i,j)); % Eq. (2.7)
                    X(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*bestPositions(j)-X(i,j)); % Eq. (2.1)
                    	X(i,j)=bestPositions(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(bestPositions(j)-X(i,j));
                % Eq. (2.5)
                X(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+bestPositions(j);
                
            end
            
        end
        
    end
end
