% The Whale Optimization Algorithm
function [Leader_pos,Convergence_curve]=WOA_DE(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=[];
FEs=0;
t=1;
% Main loop
while  FEs < MaxFEs
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        AllFitness(i) = fitness;
        FEs=FEs+1;
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,index] = sort(AllFitness,'descend');
    tmpPositions = zeros(SearchAgents_no,dim);
    for i = 1:SearchAgents_no
        tmpPositions(i,:) = Positions(index(i),:);
    end
    Positions = tmpPositions;
    re=rand();   
    reNum = floor(SearchAgents_no * re + 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%引入差分算法的交叉变异，选择机制%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = rand() * (1 - FEs / MaxFEs)^2;
    for i=1:reNum
        tmpP(i,:) = Positions(i,:);
        r = randperm(SearchAgents_no,3);
        V = Positions(r(1),:) + rand() * (Positions(r(2),:) - Positions(r(3),:));
        [OLPositions(i,:),Leader_pos,Leader_score,FEs] = DEfun(Positions(i,:),V,Leader_pos,Leader_score,f,dim,FEs,fobj);
    end
    OLP = [tmpP;OLPositions];
    OLPfitness = [];
    for i = 1:size(OLP,1)
        fitness = fobj(OLP(i,:));
        FEs = FEs + 1;
        OLPfitness = [OLPfitness fitness];
    end
    [~,index] = sort(OLPfitness);
    Positions(1:reNum,:) = OLP(index(1:reNum),:); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    a=2-FEs*((2)/MaxFEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)  
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+FEs*((-1)/MaxFEs);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                
            end
            
        end
    end
    Convergence_curve(t)=Leader_score;
    t=t+1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%差分算法%%%%%%%%%%%%%%%%%%%%%%%%
function [VPosition,Leader_pos,Leader_score,FEs] = DEfun(Position,V,Leader_pos,Leader_score,f,dim,FEs,fobj)
    CR = 0.2;
    j = randperm(dim,1);
    for ant = 1:dim
        if rand() <= CR || ant == j
            Position(ant) = V(ant);
        else
            Position(ant) = Position(ant);
        end
    end
    r = rand();
    if r >= 0.5
        VPosition = Position + f * (Position - Leader_pos);
    else
        VPosition = Position + f * (Leader_pos - Position);
    end
    Pfitness = fobj(Position);
    FEs = FEs + 1;
    if Pfitness < Leader_score
        Leader_pos = Position;
        Leader_score = Pfitness;
    end
    Vfitness = fobj(VPosition);
    FEs = FEs + 1;
    if Vfitness < Leader_score
        Leader_score = Vfitness;
        Leader_pos = VPosition;
    end
    if Vfitness < Pfitness
        VPosition = VPosition;
    else
        VPosition = Position;
    end
end


