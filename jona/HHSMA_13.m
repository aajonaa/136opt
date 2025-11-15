% Slime Mold Algorithm modified by Jona 2023-10-29.
% Update: 1.Hybrid the HHO and SMA(HHO exploration strategy applied in SMA owing to
% the better exploration performance based on the convergence_curve). 2.Use
% opposition based learning to avoid stuck in local optimal. 3.Update the
% best positition. 4.Cancel the elite opposition based learning. 5.Add the
% HHO's exploration phase to the correct position. 6.Add the HHO's
% exploitation phase 1 substitute the exploration part of HHO.
function [best_pos,convergence_curve]=HHSMA_14(N,Max_FEs,lb,ub,dim,fobj)
tic
disp('HHSMA_14 is now tackling your problem')

% initialize position
best_pos=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
convergence_curve=[];
it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

% Main loop
FEs = 0;
while  FEs < Max_FEs
    
    %sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
        FEs = FEs + 1;
    end
    
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);

    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero

    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        best_pos=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
    a = atanh(-(FEs/Max_FEs)+1);   %Eq.(2.4)
    b = 1-FEs/Max_FEs;
    
    %%%%%%%%%%
    E1=2*(1-(FEs/Max_FEs)); % factor to show the decreaing energy of rabbit
    %%%%%%%%%%
    
    % Update the Position of search agents
    for i=1:N
        % Random distribution status.
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        % Exploration status.
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(2.2) p is between [0, 1]. Decreasing from 1 to 0. 1-->0
            vb = unifrnd(-a,a,1,dim);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,dim);
            for j=1:dim
                r = rand();
                A = randi([1,N]);  % two positions randomly selected from population
                B = randi([1,N]);
                % Exploration status.
                if r<p    %Eq.(2.1)
                    X(i,j) = best_pos(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                    
                % Exploitation status.
                else
                    X(i,j) = vc(j)*X(i,j);                    
                end
            end
            
            
            % We get a new individual after the  for loop which update all
            % decision variables of the individual.
            % Then we bring the phase 2 of HHO to consider if there is a
            % between position.
            %%%%%%%%%%
            E0 = 2 * rand() - 1; %-1<E0<1
            Escaping_Energy = E1 * (E0);  % escaping energy of rabbit

%             if abs(Escaping_Energy) >= 1 % The former iteration the lager probability to be true.
%                 q = rand();
%                 rand_Hawk_index = floor(N*rand()+1);
%                 X_rand = X(rand_Hawk_index, :);
%                 if q<0.5   % perch based on other family members
%                     X1 = X_rand - rand() * abs(X_rand - 2 * rand() * X(i,:));
%                     FEs = FEs + 2;
%                     if fobj(X1) < fobj(X(i, :))
%                         X(i, :) = X1;
%                     end
% %             elseif q>=0.5   % perch on a random tall tree (random site inside group's home range)
%                 else
%                     X1 = (best_pos - mean(X)) - rand() * ((ub - lb) * rand + lb);
%                     FEs = FEs + 2;
%                     if fobj(X1) < fobj(X(i, :))
%                         X(i, :) = X1;
%                     end
%                 end
%             end
%             %%%%%%%%%%
            
            
            %%%%%%%%%%
            if abs(Escaping_Energy) < 1
                % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
                % phase 1: surprise pounce (seven kills)
                % surprise pounce (seven kills): multiple, short rapid dives by different hawks
                r = rand(); % probablity of each event
                if r >= 0.5 && abs(Escaping_Energy) < 0.5 % Hard besiege
                    X1 = (best_pos) - Escaping_Energy * abs(best_pos - X(i, :));
                    FEs = FEs + 2;
                    if fobj(X1) < fobj(X(i, :))
                        X(i, :) = X1;
                    end
                end

                if r >= 0.5 && abs(Escaping_Energy) >= 0.5  % Soft besiege
                    Jump_strength = 2 * (1 - rand()); % random jump strength of the rabbit
                    X1 = (best_pos - X(i, :)) - Escaping_Energy * abs(Jump_strength * best_pos - X(i, :));
                    FEs = FEs + 2;
                    if fobj(X1) < fobj(X(i, :))
                        X(i, :) = X1;
                    end
                end
            end
            %%%%%%%%%%
            
            
            %%%%%%%%%%
            % phase 2: performing team rapid dives (leapfrog movements)
            if abs(Escaping_Energy) < 1
            r = rand();
            if r<0.5 && abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                
                Jump_strength=2*(1-rand());
                X1=best_pos-Escaping_Energy*abs(Jump_strength*best_pos-X(i,:));
                
                FEs = FEs + 2;
                if fobj(X1)<fobj(X(i,:)) % improved move?
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=best_pos-Escaping_Energy*abs(Jump_strength*best_pos-X(i,:))+rand(1,dim).*Levy(dim);
                    
                    FEs = FEs + 2;
                    if (fobj(X2)<fobj(X(i,:))) % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1=best_pos-Escaping_Energy*abs(Jump_strength*best_pos-mean(X));
                
                FEs = FEs + 2;
                if fobj(X1)<fobj(X(i,:)) % improved move?
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit
                    X2=best_pos-Escaping_Energy*abs(Jump_strength*best_pos-mean(X))+rand(1,dim).*Levy(dim);
                    
                    FEs = FEs + 2;
                    if (fobj(X2)<fobj(X(i,:))) % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            %%%%%%%%%%
            
            
        end
    end
    
%     %%%%%%%%%%
%     dynamic_ub = max(X);
%     dynamic_lb = min(X);
%     for i = 1:N  
%         new_X(i, :) = rand() * (dynamic_ub + dynamic_lb) - X(i, :);
%         for j = 1:dim
%            if new_X(i,j) < dynamic_lb(j) || new_X(j) > dynamic_ub(j)
%                 new_X(i,j) = dynamic_lb(j) + rand() * (dynamic_ub(j) - dynamic_lb(j));  
%            end
%         end 
%     end 
%     population1 = new_X;
%     population2 = [X; population1];
%     for i = 1:2*N
%         fitness_value(i) = fobj(population2(i, :));
%         FEs = FEs + 1;
%     end
%     [~, index] = sort(fitness_value);
%     X = population2(index, :); 
%     X = X(1:N, :); 
%     %%%%%%%%%%

      % The Elite opposition based learning.
      %%%%%%%%%%
%     for i = 1:N
%         fitness_value(i) = fobj(X(i, :));
%         FEs = FEs + 1;
%     end
%     [~, index] = sort(fitness_value);
%     best_pos = X(index(1), :);
%     Temp = ub + lb - best_pos;
%     Flag4ub=Temp>ub;
%     Flag4lb=Temp<lb;
%     Temp=(Temp.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
%     fTemp = fobj(Temp);
%     FEs = FEs + 1;
%     if fTemp<bestFitness
%         bestFitness = fTemp;
%         best_pos = Temp;
%     end
      %%%%%%%%%%
    
    convergence_curve(it)=Destination_fitness;
    it=it+1;
end
toc
end

