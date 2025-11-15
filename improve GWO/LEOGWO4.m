
%  Developed in MATLAB R2013a                                    %

% Levy opposition-based learning Grey Wolf Optimizer   Levy+EOL+GWO+spiral
% LF-based motions
% Authors and programmers
% Ali Asghar Heidari-  University of Tehran- as_heidari@ut.ac.ir
% Huiling Chen,Jie Luo -University of WenZhou.

%note: update in 20170719 by Jie Luo%
%1:fixed the EOL formula  in opposition based part
%2:to boost exploration,adding movements toward the random leader and levy
%   flight operator when A>1
%
function [Alpha_pos,Convergence_curve]=LEOGWO4(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
% display('LEOGWO3 is optimizing your problem');

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

n=SearchAgents_no;
%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
Oppositions=zeros(SearchAgents_no,dim); %opposition
fitness=zeros(SearchAgents_no,1);
oppofitness=zeros(SearchAgents_no,1);

%Convergence_curve=zeros(1,Max_iter);
Convergence_curve=[];
Fes=0;% Loop counter
l = 1;
% Main loop
while Fes<Max_iter
    %
    for i=1:size(Positions,1)
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
		if Fes < Max_iter
			Fes = Fes + 1;
			% Calculate objective function for each search agent
			fitness(i)=fobj(Positions(i,:));
			betterfitness=fitness(i);
			betterPos=Positions(i,:);
			
			% Update Alpha, Beta, and Delta
			if betterfitness<Alpha_score
				Alpha_score=betterfitness; % Update alpha
				Alpha_pos=betterPos;
			end
			if betterfitness>Alpha_score && betterfitness<Beta_score
				Beta_score=betterfitness; % Update beta
				Beta_pos=betterPos;
			end
        else
			break;
		end
    end
    %% opposition based part

    for i=1:size(Positions,1),
        Oppositions(i,:)=lb+ub-Alpha_pos+rand.*(Alpha_pos-Positions(i,:));
        % Calculate objective function for each search agent
        Flag4ub=Oppositions(i,:)>ub;
        Flag4lb=Oppositions(i,:)<lb;
        Oppositions(i,:)=(Oppositions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
		if Fes < Max_iter
			Fes = Fes + 1;
			oppofitness(i)=fobj(Oppositions(i,:));
			
			betterfitness=oppofitness(i);
			betterPos=Oppositions(i,:);
			
			% Update Alpha, Beta, and Delta
			if betterfitness<Alpha_score
				Alpha_score=betterfitness; % Update alpha
				Alpha_pos=betterPos;
			end
			
			if betterfitness>Alpha_score && betterfitness<Beta_score
				Beta_score=betterfitness; % Update beta
				Beta_pos=betterPos;
			end
        else
			break;
		end
    end
    
    AllPositions=[Positions;Oppositions];
    AllFitness=[fitness;oppofitness];
    % sort gwo
    [~,Index]=sort(AllFitness);
    %     [AllFitness,Index]=sort(AllFitness,'descend');
    AllPositions_tmp=AllPositions;
    for i=1:size(AllPositions,1),
        AllPositions(i,:)=AllPositions_tmp(Index(i),:);
    end
    
    Positions=AllPositions(1:n,:);%将原来灰狼和反向学习后的灰狼按适应度值排序后的前n个灰狼的位置
    %% the rest
    
    a=2-Fes*((2)/Max_iter); % a decreases linearly fron 2 to 0
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+Fes*((-1)/Max_iter);
    for i=1:size(Positions,1)
        A=2*a*rand-a;  % Eq. (2.3) in the paper
        A1=2*a*rand(1,dim)-a; % Equation (3.3)
        C1=2*rand(1,dim); % Equation (3.4)
        
        A2=2*a*rand(1,dim)-a; % Equation (3.3)
        C2=2*rand(1,dim); % Equation (3.4)

        b=1;               %  parameters in Eq. (2.5)
        t=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        if p<0.5
            beta=1*rand()*a;
            D_alpha=abs(C1.*Alpha_pos-Positions(i,:));
            D_beta=abs(C2.*Beta_pos-Positions(i,:)); % Equation (3.5)-part 2
            if abs(A)>=1
                 %move to the random leader
                rand_leader_index = floor(SearchAgents_no*rand()+1);
                X_rand = Positions(rand_leader_index, :);
                D_X_rand=abs(C1.*X_rand-Positions(i,:)); % Eq. (2.7)
                newPositions1=X_rand-A1.*D_X_rand;  
                 % Return back the search agents that go beyond the boundaries of the search space
                Flag4ub=newPositions1>ub;
                Flag4lb=newPositions1<lb;
                newPositions1=(newPositions1.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                % levy flight
                newPositions2=0.5*((Alpha_pos-A1.*D_alpha)+...
                                 (Beta_pos-A2.*D_beta))+...
                                  rand(1,dim).*Levy(dim,beta);
                % Return back the search agents that go beyond the boundaries of the search space
                Flag4ub=newPositions2>ub;
                Flag4lb=newPositions2<lb;
                newPositions2=(newPositions2.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                %GS
                if fobj(newPositions1)<fobj(newPositions2)
                    Positions(i,:)=newPositions1;
                else
                    Positions(i,:)=newPositions2;
                end
            elseif abs(A)<1
                    % move to leader
                     newPositions=0.5*((Alpha_pos-A1.*D_alpha)+...
                                          (Beta_pos-A2.*D_beta));
                    % Return back the search agents that go beyond the boundaries of the search space
                    Flag4ub=newPositions>ub;
                    Flag4lb=newPositions<lb;
                    newPositions=(newPositions.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                    %GS
                    if fobj(newPositions)<fobj(Positions(i,:))
                        Positions(i,:)=newPositions;
                    end
            end            
    %% ------ spiral levy-embedded part: it move based on a new rule like WOA-------------        
        elseif p>=0.5
            %random beta
             beta=1*rand()*a;
            
            distance2Leader1=abs(Alpha_pos-Positions(i,:));
            distance2Leader2=abs(Beta_pos-Positions(i,:));
            
            % Eq. (2.5)
            Positions(i,:)=1*(Alpha_pos)...
            + distance2Leader1*exp(b.*t).*cos(t.*2*pi) ...
            + distance2Leader2*exp(b.*t).*cos(t.*2*pi)...    
            + rand(1,dim).*Levy(dim,beta);
        end
    %% ------ end of spiral levy-embedded part -------------        
    end
    Convergence_curve(l)=Alpha_score;
    l=l+1;
    
end
end

function o=Levy(d,beta)

%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);

% Eq. (3.9)
o=step;

end


