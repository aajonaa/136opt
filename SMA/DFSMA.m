% Slime Mold Algorithem
function [bestPositions,Convergence_curve]=DFSMA(N,MaxFEs,lb,ub,dim,fobj)

% initialize position
bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
%X1=initialization(N,dim,ub,lb);
Convergence_curve=[];
it=1;
FEs=0;
X1 = X;%%%%%%%%%%%
z = 0.03;%%%%%%%%%%%
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

S=bestFitness-worstFitness+eps;%%%%%%%%%%%%
	
if bestFitness < Destination_fitness
        bestPositions=X1(SmellIndex(1),:);
        Destination_fitness = bestFitness;
end
    
% Main loop
while  FEs < MaxFEs    
    %calculate the fitness weight of each slime mold
    for i=1:size(X,1)
        for j=2:size(X,2)
            if i<=(N/2)    %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
            
        end
    end
    
    a = atanh(-(FEs/MaxFEs)+1); %Eq.(2.4)%%%%%%%%%%%
    b = 1-FEs/MaxFEs; % a decreases linearly fron 2 to 0	%%%%%%%%%%%
	pa=1-(1-0.6)*FEs/MaxFEs;%%%%%%%%%%%
    freq=0.4-(0.4-0)*FEs/MaxFEs;	%%%%%%%%%%%
	
    % Update the Position of search agents
    for i=1:N
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        else%%%%%%%%%%%
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
    
   for i=1:N
		Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
		FEs = FEs + 1;
		temp_fit = fobj(X(i,:));
	    
        if AllFitness(i) > temp_fit
			AllFitness(i) = temp_fit;
			X1(i,:) = X(i,:);
        end
   end   
	[SmellOrder,SmellIndex] = sort(AllFitness);
    worstFitness = SmellOrder(size(X,1));
    bestFitness = SmellOrder(1);
    S=bestFitness-worstFitness+eps;
    
	if bestFitness < Destination_fitness
        bestPositions=X1(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
	nest=X1;
    K=rand(N,dim)>freq;
    stepsize=normrnd(0.1,0.5)*(nest(randperm(N),:)-nest(randperm(N),:));
    temp_X =nest+stepsize.*K;
	for i=1:N
		Flag4ub=temp_X(i,:)>ub;
        Flag4lb=temp_X(i,:)<lb;
        temp_X(i,:)=(temp_X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
		FEs = FEs + 1;
		temp_fit = fobj(temp_X(i,:));
		if AllFitness(i) > temp_fit
			AllFitness(i) = temp_fit;
			X1(i,:) = temp_X(i,:);
        end	
    end
        [SmellOrder,SmellIndex] = sort(AllFitness);
        worstFitness = SmellOrder(size(X,1));
        bestFitness = SmellOrder(1);
        S=bestFitness-worstFitness+eps;
	
        if bestFitness < Destination_fitness
            bestPositions=X1(SmellIndex(1),:);
            Destination_fitness = bestFitness;
        end     
        
        
        Convergence_curve(it)=Destination_fitness;
        %display();
        it=it+1;
end
end



