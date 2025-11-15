% The Whale Optimization Algorithm
%% Gauss+Cauzy+混沌初始化+CLS 2018.3.1 by CHL
function [Leader_pos,Convergence_curve]=CCMWOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems
K = SearchAgents_no;

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
m = 1500;
FEs=0;

%%混沌初始
x = rand();
if(x~=0.25&&x~=0.5&&x~=0.75)
    ch(1) = x;
end
for ant=1:(SearchAgents_no-1)
    ch(ant+1)=4*ch(ant)*(1-ch(ant));
    PCh = ch(ant)*Positions; 
    PHe = [Positions;PCh];
    count=size(PHe,1);
    FitnessHe1=[];
    for i=1:count
        PHeLin=fobj(PHe(i,:));
        FEs=FEs+1;
        FitnessHe1 = [FitnessHe1 PHeLin];
    end
%     FitnessHe1=calculateFitness(ObjValHe1);
    [FitnessHe2,index] = sort(FitnessHe1);
    Positions = PHe(index,:);
    Positions = Positions(1:SearchAgents_no,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Convergence_curve=[];
% FEs=0;
t=1;
% Main loop
while  FEs < MaxFEs
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        if FEs<MaxFEs
            FEs=FEs+1;
            fitness=fobj(Positions(i,:));
            % Update the leader
            if fitness<Leader_score % Change this to > for maximization problem
                Leader_score=fitness; % Update alpha
                Leader_pos=Positions(i,:);
            end
        else 
            break;
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
            
            Moth_pos_m_gaus=Positions(i,:)*(1+randn(1));
            if FEs+1<MaxFEs
                FEs=FEs+2;
                Moth_fitness_m_gaus=fobj(Moth_pos_m_gaus);
                Moth_fitness_s=fobj(Positions(i,:));
                if Moth_fitness_m_gaus<Moth_fitness_s        
                    Positions(i,:)=Moth_pos_m_gaus; 
                end
            else 
                break;
            end
        end
   end
    
        setCan = 1-power(abs((FEs-1)/FEs),m);

%     setCan = (MaxFEs-t+1)/MaxFEs;
    x = rand();
    if(x~=0.25&&x~=0.5&&x~=0.75)
        ch(1) = x;
    end
    for ant=1:(K)
        ch(ant+1)=4*ch(ant)*(1-ch(ant));
        CH(ant,:) = lb+ch(ant)*(ub-lb);    %ub大
        V = (1-setCan)*Leader_pos+setCan*CH(ant);
        if FEs<MaxFEs
            FEs=FEs+1;
            [FitnessV]=fobj(V);%计算适应度值
            if (FitnessV<Leader_score)
                Leader_score = FitnessV;
                Leader_pos = V;
                break;
            end
        else
            break;
        end
    end
    Convergence_curve(t)=Leader_score;
    t=t+1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




