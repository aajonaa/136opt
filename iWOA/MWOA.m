
%% Sun, Y., et al. (2018). "A modified whale optimization algorithm for large-scale global optimization problems." Expert Systems with Applications 114: 563-577.
%10.1016/j.eswa.2018.08.027
%% WOA+QI+Levy
function [Leader_pos,Convergence_curve]=MWOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

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
        fitness(i,:)=fobj(Positions(i,:));
        FEs=FEs+1;
        % Update the leader
        if fitness(i,:)<Leader_score % Change this to > for maximization problem
            Leader_score=fitness(i,:); % Update alpha
            Leader_pos=Positions(i,:);
            Leader_index=i;
        end
        
    end
    
    a=2*cos(FEs/MaxFEs); %nonlinear control parameter strategy
%     a=2-FEs*((2)/MaxFEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
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
        p2=rand();
        
        k=randperm(SearchAgents_no,2);
        while k(1)==Leader_index || k(2)==Leader_index
            k=randperm(SearchAgents_no,2);
        end
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    %% Levy 
                    beta=unifrnd(0,2);
                    P(i,j)=Levy(beta,Positions(i,j),Leader_pos(j));
                    Positions(i,j)=Positions(i,j)+1/sqrt(t)*sign(rand()-0.5)*P(i,j);
%                     D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
%                     Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
                if p2<0.6
                    distance2Leader=abs(Leader_pos(j)-Positions(i,j));
%                     Eq. (2.5)
                    Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                    
                elseif p2>=0.6
                    %% QI
                    member=(Positions(k(1),j)^2-Positions(k(2),j)^2)*fitness(Leader_index)+(Positions(k(2),j)^2-Leader_pos(j)^2)*fitness(k(1))+(Leader_pos(j)^2-Positions(k(1),j)^2)*fitness(k(1));
                    denominator=(Positions(k(1),j)-Positions(k(2),j))*fitness(Leader_index)+(Positions(k(2),j)-Leader_pos(j))*fitness(k(1))+(Leader_pos(j)-Positions(k(1),j))*fitness(k(1));
                    Positions(i,j)=0.5*member/denominator;
                end    
                
            end
            
        end
    end
    Convergence_curve(t)=Leader_score;
    t=t+1;

end

end
function P=Levy(beta,X,Xbest) %% d是[0,2]之间的随机数
%     beta=d;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn()*sigma;
    v=randn();
    step=u/abs(v)^(1/beta);
    L=0.01*step*(X-Xbest);
    P=randn()*L;
end


