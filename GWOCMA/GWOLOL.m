
%用beta来初始化种群，前期在levy和GWO中贪心选择最优，后期用CMA
% Grey Wolf Optimizer
function [Alpha_pos,Convergence_curve]=GWOLOL(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems
Q=3;
F=4;

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
N = SearchAgents_no;
it=1;FEs=0;
Convergence_curve=[];
R = 0.5;
flag = 1;
for i=1:SearchAgents_no
    Positions(i,:) = Positions(i,:) *(1+betarnd(1.2,1.2));
end
%cov = covariance(Positions);
% Main loop
while  FEs < MaxFEs
    % % Main loop
    % while l<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        tag = 1;
        if FEs<MaxFEs
            FEs=FEs+1;
            % Calculate objective function for each search agent
            fitness(1,i)=fobj(Positions(i,:));
            
            % Update Alpha, Beta, and Delta
            if fitness(1,i)<Alpha_score
                Alpha_score=fitness(1,i); % Update alpha
                Alpha_pos=Positions(i,:);
                alpha_no = i;
                if tag == 1
                    beta_no1 = alpha_no;
                    delta_no = alpha_no;
                    tag = 0;
                end
            end
            
            if fitness(1,i)>Alpha_score && fitness(1,i)<Beta_score
                Beta_score=fitness(1,i); % Update beta
                Beta_pos=Positions(i,:);
                beta_no1 = i;
            end
            
            if fitness(1,i)>Alpha_score && fitness(1,i)>Beta_score && fitness(1,i)<Delta_score
                Delta_score=fitness(1,i); % Update delta
                Delta_pos=Positions(i,:);
                delta_no = i;
            end
        else
            break;
        end
    end
    % Update the Position of search agents including omegas           
        a=2-FEs*((2)/MaxFEs); % a decreases linearly fron 2 to 0
        for i=1:size(Positions,1)        
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
%                 w1 = Alpha_score/(Alpha_score + Beta_score + Delta_score);
%                 w2 = Beta_score/(Alpha_score + Beta_score + Delta_score);
%                 w3 = Delta_score/(Alpha_score + Beta_score + Delta_score);
%                 Positions(i,j) = X1 * w1 + X2 * w2 + w3 * X3; 
%                   temp(1,j) =(X1*w1+X2*w2+X3*w3)/3;% Equation (3.7)
%                 temp(1,j)=(X1+X2+X3)/3;% Equation (3.7) 
              Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)           
         end
%            if fobj(temp) <= fobj(temp1)
%                Positions(i,:) = temp;
%            else
%                Positions(i,:) = temp1;
%            end

		[ k1,k2,k3 ] = GetRan3( i,N );
		T=Positions(k1,:)+rand(1,dim).*(Positions(k2,:)-Positions(k3,:));		
        [ OLPosition(i,:),OLFitness(1,i),FEs ] = Con_OEDCanSolution2( Positions(i,:),T,dim,Q,F,fobj,FEs );
		if OLFitness(1,i)<fitness(1,i)
			Positions(i,:)=OLPosition(i,:);
			fitness(1,i)=OLFitness(1,i);
        end
        
        end
		X3=Lenvy(Positions,1);
        [Positions,FEs]=Best_Position(Positions,X3,fobj,FEs); 
    
     
        
        %     Alpha_score = sorted_fit(1);
        %     Alpha_pos = Positions(1,:);
    
    %     FEs = FEs + 1;
    %     Convergence_curve(FEs)=Alpha_score;
    Convergence_curve(it)=Alpha_score;
    it=it+1;
end
end
function o=L(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end


function [X,FEs]=Best_Position(X,X3,fobj,FEs)
n=size(X,1);
for i=1:n
	FEs = FEs + 2;
    y_X=fobj(X(i,:)); 
    y_X3=fobj(X3(i,:));
    y_col=[y_X,y_X3];
    [~,k]=min(y_col);
    switch k
        case 1
            X(i,:)=X(i,:);
        case 2
            X(i,:)=X3(i,:);       
    end
end               
end




function X=Lenvy(X0,k)
n=size(X0,1);
for i=1:n
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn*sigma;
    v=randn;
    step=u./abs(v).^(1/beta);
    o=step;
    X1(i,:)=X0(i,:)*(1+k*o);
end
X=X1;
end