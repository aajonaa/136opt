%% Wu, J., et al. (2018). "Path planning for solar-powered UAV in urban environment." Neurocomputing 275: 2055-2065.
%10.1016/j.neucom.2017.10.037
%% WOA+adaptive chaos-gaussian switching solving strategy+coordinated decision-making strategy
function [Leader_pos,Convergence_curve]=IWOA2(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems
FEs=0;
%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);


Convergence_curve=[];
Nb=0;
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
    
    %% rank the fitness,obtain the top three optimal values
    [~,index]=sort(fitness);
    X1=Positions(index(1),:);
    X2=Positions(index(2),:);
    X3=Positions(index(3),:);
    
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
            
            %% Coordinated decision-making mechanism
            if p<0.5   
                D=abs(C*(X1(j)+X2(j)+X3(j))/3-Positions(i,j));
                Positions(i,j)=(X1(j)+X2(j)+X3(j))/3-A*D;
%                 if abs(A)>=1
%                     rand_leader_index = floor(SearchAgents_no*rand()+1);
%                     X_rand = Positions(rand_leader_index, :);
%                     D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
%                     Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
%                     
%                 elseif abs(A)<1
%                     D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
%                     Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
%                 end
                
            elseif p>=0.5
              
%                 distance2Leader=abs(Leader_pos(j)-Positions(i,j));
%                 % Eq. (2.5)
%                 Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                Dp=abs((X1(j)+X2(j)+X3(j))/3-Positions(i,j));
                Positions(i,j)=(X1(j)+X2(j)+X3(j))/3+Dp*exp(b.*l).*cos(l.*2*pi);
            end
            
        end
        newfitness(i,:)=fobj(Positions(i,:));
        FEs=FEs+1;
        if newfitness(i,:)<fitness(i,:)
            Nb=Nb+1;
        end
    end
    
    %% ADAPTIVE CHAOS-GAUSSIAN SWITCHING SLOVINE STRATEGY
    IR=Nb/SearchAgents_no;
    xmax=max(Positions);
    xmin=min(Positions);
    beta=rand(); %chaos parameter
    u=4;
    zita=1-abs((FEs-1)/FEs);
    while beta==0.25 || beta==0.5 || beta==0.75
        beta=rand();
    end
    for i=1:size(Positions,1)
        for j=1:dim
            % Gaussian variation
            if IR<0.15
                delta_g=xmin(j)+randn()*(xmax(j)-xmin(j));
                Positions(i,j)=(1-zita)*Positions(i,j)+zita*delta_g;
            end
            
            %no variation
            if IR>=0.15 && IR<=0.25
                Positions(i,j)=Positions(i,j);            
            end
            
            % chaos variation
            if IR>0.25
                beta=u*beta*(1-beta);
                delta_c=xmin(j)+beta*(xmax(j)-xmin(j));
                Positions(i,j)=(1-zita)*Positions(i,j)+zita*delta_c;
            end
        end
    end

        
    Convergence_curve(t)=Leader_score;
    t=t+1;

end

