%% The Whale Optimization Algorithm
%% Abd Elaziz, M. and D. Oliva (2018). "Parameter estimation of solar cells diode models by an improved opposition-based whale optimization algorithm." Energy Conversion and Management 171: 1843-1859.
% 10.1016/j.enconman.2018.05.062 
function [Leader_pos,Convergence_curve]=OBWOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
FEs=0;
for i=1:size(Positions,1)
    fitness(i,:)=fobj(Positions(i,:));
    FEs=FEs+1;
end

%% OBL-initializal agents 一般反向学习 静态边界
[P,fes,newfitness]=OBL(Positions,fitness,SearchAgents_no,dim,ub,lb,fobj,1);
Positions=P;
FEs=FEs+fes;
fitness=newfitness;
Leader_score=min(fitness);

Convergence_curve=[];

t=1;
% Main loop
while  FEs < MaxFEs
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        fitness(i,:)=fobj(Positions(i,:));
        FEs=FEs+1;
        % Update the leader
        if fitness(i,:)<Leader_score % Change this to > for maximization problem
            Leader_score=fitness(i,:); % Update alpha
            Leader_pos=Positions(i,:);
        end
        
    end
    
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
%% OBL--一般反向学习 静态边界
[P,fes,newfitness]=OBL(Positions,fitness,SearchAgents_no,dim,ub,lb,fobj,1);
Positions=P;
FEs=FEs+fes;
fitness=newfitness;
Leader_score=min(fitness);

    Convergence_curve(t)=Leader_score;
    t=t+1;

end
end





