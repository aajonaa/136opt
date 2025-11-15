

% Grey Wolf Optimizer
function [Alpha_pos,Convergence_curve]=CAGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems
Fitness=zeros(1,SearchAgents_no);
Temp=zeros(1,dim);
type=2;
%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=[];

Fes=0;% Loop counter
it = 1;
% Main loop
while Fes<Max_iter
     for i=1:size(Positions,1)  
%         
%        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
%         
%         % Calculate objective function for each search agent
%         fitness=fobj(Positions(i,:));
%         
%         % Update Alpha, Beta, and Delta
%         if fitness<Alpha_score 
%             Alpha_score=fitness; % Update alpha
%             Alpha_pos=Positions(i,:);
%         end
%         
%         if fitness>Alpha_score && fitness<Beta_score 
%             Beta_score=fitness; % Update beta
%             Beta_pos=Positions(i,:);
%         end
%         
%         if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
%             Delta_score=fitness; % Update delta
%             Delta_pos=Positions(i,:);
%         end
     end
    
    
    a=2-Fes*((2)/Max_iter); % a decreases linearly fron 2 to 0
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
		Fes = Fes + 1;
        Fitness(1,i)=fobj(Positions(i,:));
        [ Alpha_pos,Beta_pos,Delta_pos ,Fes ] = GCA( i,Positions,SearchAgents_no,dim,fobj,type,Fes ); %i：行数（第几个灰狼）Positions：灰狼位置信息 SearchAgents_no：灰狼个数 dim：维度 fobj：目标函数  type：类型
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Temp(j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
		Fes = Fes + 1;
        TempFitness=fobj(Temp);
        if TempFitness<Fitness(1,i)
            Fitness(1,i)=TempFitness;
            Positions(i,:)=Temp;
        end
        if Fitness(1,i)<Alpha_score
            Alpha_score=Fitness(1,i);
        end
    end
    Convergence_curve(it)=Alpha_score;
	it=it+1;   
end



