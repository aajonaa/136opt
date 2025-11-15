% % The Whale Optimization Algorithm
% function [best_pos,Convergence_curve]=WOA(N,Max_FEs,lb,ub,dim,fobj)
% 
%     % initialize position vector and score for the leader
%     best_pos=zeros(1,dim);
%     bestFitness=inf; %change this to -inf for maximization problems
%     AllFitness = inf * ones(N, 1);
% 
% 
%     %Initialize the positions of search agents
%     X=initialization(N,dim,ub,lb);
% 
%     Convergence_curve=[];
%     FEs=0;
%     t=1;
%     % Main loop
%     while  FEs < Max_FEs
%         for i=1:N
% 
%             % Return back the search agents that go beyond the boundaries of the search space
%             Flag4ub=X(i,:)>ub;
%             Flag4lb=X(i,:)<lb;
%             X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
% 
%             % Calculate objective function for each search agent
% %             fitness=fobj(X(i,:));
%             AllFitness(i, 1)=fobj(X(i,:));
%             FEs=FEs+1;
%             % Update the leader
% %             if fitness<bestFitness % Change this to > for maximization problem
%             if AllFitness(i, 1)<bestFitness % Change this to > for maximization problem
% %                 bestFitness=fitness; % Update alpha
%                 bestFitness=AllFitness(i, 1); % Update alpha
%                 best_pos=X(i,:);
%             end
% 
%         end
%         
%         [sorted_AllFitness, ~] = sort(AllFitness);
% 
%         a=2-FEs*((2)/Max_FEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)
% 
%         % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
%         a2=-1+FEs*((-1)/Max_FEs);
% 
%         Rolette_index=RouletteWheelSelection(1./sorted_AllFitness);
%         if Rolette_index==-1  
%             Rolette_index=1;
%         end
% 
%         % Update the Position of search agents 
%         for i=1:N
%             r1=rand(); % r1 is a random number in [0,1]
%             r2=rand(); % r2 is a random number in [0,1]
% 
%             A=2*a*r1-a;  % Eq. (2.3) in the paper
%             C=2*r2;      % Eq. (2.4) in the paper
% 
% 
%             b=1;               %  parameters in Eq. (2.5)
%             l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
% 
%             p = rand();        % p in Eq. (2.6)
% 
%             for j=dim
% 
%                 if p<0.5   
%                     if abs(A)>=1
%                         rand_leader_index = floor(N*rand()+1);
%                         X_rand = X(rand_leader_index, :);
%                         D_X_rand=abs(C*X_rand(j)-X(i,j)); % Eq. (2.7)
%                         X(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
% 
%                     elseif abs(A)<1
%                         D_Leader=abs(C*X(Rolette_index, j)-X(i,j)); % Eq. (2.1)
%                         X(i,j)=best_pos(j)-A*D_Leader;      % Eq. (2.2)
%                     end
% 
%                 elseif p>=0.5
%     %               
%                     distance2Leader=abs(best_pos(j)-X(i,j));
%                     % Eq. (2.5)
%                     X(i,j)=best_pos(j) + distance2Leader*exp(b.*l).*cos(l.*2*pi);
% 
%                 end
% 
%             end
%         end
%         Convergence_curve(t)=bestFitness;
%         t=t+1;
% 
%     end
% end

% The Whale Optimization Algorithm
function [best_pos,Convergence_curve]=WOA(N,Max_FEs,lb,ub,dim,fobj)
tic

% initialize position vector and score for the leader
best_pos=zeros(1,dim);
bestFitness=inf; %change this to -inf for maximization problems
AllFitness = inf * ones(N, 1);


%Initialize the positions of search agents
X=initialization(N,dim,ub,lb);

Convergence_curve=[];
FEs=0;
t=1;
% Main loop
while  FEs < Max_FEs
    for i=1:N
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
%         fitness=fobj(X(i,:));
        AllFitness(i, 1)=fobj(X(i,:));
        FEs=FEs+1;
%             if fitness<bestFitness % Change this to > for maximization problem
        if AllFitness(i, 1)<bestFitness % Change this to > for maximization problem
%                 bestFitness=fitness; % Update alpha
            bestFitness=AllFitness(i, 1); % Update alpha
            best_pos=X(i,:);
        end
        
    end
    
    [sorted_AllFitness, ~] = sort(AllFitness);
    
    a=2-FEs*((2)/Max_FEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+FEs*((-1)/Max_FEs);
    
    Rolette_index=RouletteWheelSelection(1./sorted_AllFitness);
    if Rolette_index==-1  
        Rolette_index=1;
    end
    
    % Update the Position of search agents 
    for i=1:N
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(X,2)
            
%             if p<0.5
            if rand < rand  
                if abs(A)>=1 % global search
%                     rand_leader_index = floor(N*rand()+1);
%                     X_rand = X(rand_leader_index, :);
%                     X_rand = X(Rolette_index, :);
%                     D_X_rand=abs(C*X_rand(j)-X(i,j)); % Eq. (2.7)
%                     X(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    X(i, j) = X(i, j) + A * randn * (X(i, j) - X(floor(N*rand()+1), j));
                    
                elseif abs(A)<1 % exploitation
%                     D_Leader=abs(C*best_pos(j)-X(i,j)); % Eq. (2.1)
%                     D_Leader=abs(C*X(Rolette_index, j)-X(i,j)); % Eq. (2.1)
%                     X(i,j)=best_pos(j)-A*D_Leader;      % Eq. (2.2)negn
%                     X(i,j)=X(Rolette_index, j)-A*D_Leader;      % Eq. (2.2)
                    X(i,j) = X(Rolette_index, j) + A * abs(X(i,j) - C * X(Rolette_index, j));      % Eq. (2.2)
                end
                
%             elseif p>=0.5
            else
              
%                 distance2Leader=abs(best_pos(j)-X(i,j));
                % Eq. (2.5)
                X(i,j)=best_pos(j) + abs(X(i,j) - best_pos(j)) * exp(b*l) * cos(l*2*pi); 
                
            end
            
        end
    end
    Convergence_curve(t)=bestFitness;
    t=t+1;

end
toc
end

function choice = RouletteWheelSelection(weights)
  accumulation = cumsum(weights);
  p = rand() * accumulation(end);
  chosen_index = -1;
  for index = 1:length(accumulation)
    if (accumulation(index) > p)
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
end



