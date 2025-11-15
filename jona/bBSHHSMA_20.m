% Slime Mold Algorithm binary version modified by Jona 2023-10-29.
% Update: 1.Hybrid the HHO and SMA(HHO exploration strategy applied in SMA owing to
% the better exploration performance based on the convergence_curve). 2.Use
% opposition based learning to avoid stuck in local optimal. 3.Update the
% best positition. 4.Cancel the elite opposition based learning. 5.Add the
% HHO's exploration phase to the correct position. 6.Add the HHO's
% exploitation phase 1 substitute the exploration part of HHO. 7.Add the
% HHO's exploitation phase2. 8.Cacel 7 and add the BSA for worst 5
% individuals.9.Only use BSA in the rear part of the evalution test. 10.
% The worst 10 individual applied to BSA. 11.Change the search agents size
% to 15. 12.Change the search agents size to 30 at 1/3 FEs and canceled the
% HHO.
function [Destination_fitness, best_pos,convergence_curve, Time]=bBSHHSMA_20(N,Max_iter,dim,A, trn, vald, classifierFhd)
tic
disp('bBSHHSMA_20 is now tackling your problem')

TFid = 8;

% initialize position
best_pos=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
% X=initialization(N,dim,ub,lb); 
X=initialization(N,dim,ub,lb) > 0.5;
convergence_curve=[];
% it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

%%%%%%%%%%BSA
% bs_X = X([1, 2, 3, 4, 5], :);
for i = 1:30
    bs_X(i, :) = X(i, :);
end
bs_h_X = X;
%%%%%%%%%%

% Main loop
FEs = 0;
while  FEs < Max_iter
    
    %sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
%         Flag4ub=X(i,:)>ub;
%         Flag4lb=X(i,:)<lb;
%         X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = AccSz2(X(i, :), A, trn, vald, classifierFhd);
%         FEs = FEs + 1;
    end

    [sorted_fitness_value, sorted_fitness_index] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = sorted_fitness_value(N);
    bestFitness = sorted_fitness_value(1);
    
    %%%%%%%%%%
    if FEs < Max_iter / 3
        % The worst five individuals and its fitness.
%         bs_X_fitness = inf * ones(1, 5);
        bs_X_fitness = inf * ones(1, 30);
%         bs_X_fitness = sorted_fitness_value([N - 4, N  - 3, N - 2, N - 1, N]);
%         bs_X = X([sorted_fitness_index(N-4), sorted_fitness_index(N-3), sorted_fitness_index(N-2), sorted_fitness_index(N-1), sorted_fitness_index(N)], :);
        j = 1;
        for i = 1:N
            bs_X_fitness(j) = sorted_fitness_value(i, :);
            bs_X(j, :) = X(sorted_fitness_index(i), :);
            j = j + 1;
        end
%         bs_h_X = bs_X(randperm(5), :);
        bs_h_X = bs_X(randperm(N), :);
        F = 3 * randn;
%         map = zeros(5, dim);
        map = zeros(N, dim);
%         for i = 1:5
%             u = randperm(dim);
%             map(i, u(1:ceil(rand * dim))) = 1;
%         end
        for i = 1:30
            u = randperm(dim);
            map(i, u(1:ceil(rand * dim))) = 1;
        end

    %     for i = 1:5
    %         map(i, randi(dim)) = 1;
    %     end

        offsprings = bs_X + (map .* F) .* (bs_h_X - bs_X);
%         offsprings = BoundaryControl(offsprings, lb, ub);
        for i = 1:N
            for j = 1:dim
                offsprings(i, j) = transferFun(offsprings(i, j), offsprings(i, j), TFid);
            end
        end

        % Update the Destination fitness and the best individual after generate
        % new generation.
%         fitness_offsprings = inf * ones(1, 5);
        fitness_offsprings = inf * ones(1, N);
%         for i = 1:5
%             fitness_offsprings(i) = AccSz2(offsprings(i, :), A, trn, vald, classifierFhd);
%             FEs = FEs + 1;
%             if fitness_offsprings(i) < bestFitness
%                 bestFitness = fitness_offsprings(i);
%                 best_pos = offsprings(i, :);
%             end
%         end    

        %%%% Useless code?
        for i = 1:N
            fitness_offsprings(i) = AccSz2(offsprings(i, :), A, trn, vald, classifierFhd);
%             FEs = FEs + 1;
%             if fitness_offsprings(i) < bestFitness
%                 bestFitness = fitness_offsprings(i);
%                 best_pos = offsprings(i, :);
%             end
        end
        %%%%

        % To choose the best 5 form two generation.
        bs_X_fitness = reshape(bs_X_fitness, [1, N]);
        ind = fitness_offsprings < bs_X_fitness;
    %     disp('fitness_offsprings');
    %     disp(fitness_offsprings);
    %     disp('bs_X_fitness');
    %     disp(bs_X_fitness);
    %     disp('ind');
    %     disp(ind);
        bs_X_fitness(ind) = fitness_offsprings(ind);
        bs_X(ind, :) = offsprings(ind, :);

        %%%%Useless
    %     [minimum_fitness, index] = min(bs_X_fitness);
    %     maximum_fitness = max(bs_X_fitness);
    %     if minimun_fitness < bestFitness
    %         bestFitness = minimum_fitness;
    %         best_pos = bs_X(index, :);
    %     end
    %     if maximum_fitness > worstFitness
    %         worstFitness = maximum_fitness;
    %     end
        %%%%

%         X([N-4, N-3, N-2, N-1, N], :) = bs_X;
%         AllFitness([N-4, N-3, N-2, N-1, N]) = bs_X_fitness;
        j = 1;
        for i = 1:N
            X(i, :) = bs_X(j, :);
            AllFitness(i) = bs_X_fitness(j);
            j = j + 1;
        end

        % Resort the fitness of X.
        [sorted_fitness_value, sorted_fitness_index] = sort(AllFitness);  %Eq.(2.6)
        worstFitness = sorted_fitness_value(N);
        bestFitness = sorted_fitness_value(1);
    end
    %%%%%%%%%%

    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero

    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(2.5)
                weight(sorted_fitness_index(i),j) = 1+rand()*log10((bestFitness-sorted_fitness_value(i))/(S)+1);
            else
                weight(sorted_fitness_index(i),j) = 1-rand()*log10((bestFitness-sorted_fitness_value(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        best_pos=X(sorted_fitness_index(1),:);
        Destination_fitness = bestFitness;
    end
    
    a = atanh(-(FEs/Max_iter)+1);   %Eq.(2.4)
    b = 1-FEs/Max_iter;
    
%     %%%%%%%%%%HHO
%     E1=2*(1-(FEs/Max_iter)); % factor to show the decreaing energy of rabbit
%     %%%%%%%%%%
    
    % Update the Position of search agents
    for i=1:N
        % Random distribution status.
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
            for j = 1:dim
                X(i, j) = transferFun(X(i, j), X(i, j), TFid);
            end
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
                    X(i, j) = transferFun(X(i, j), X(i, j), TFid);
                    
                % Exploitation status.
                else
                    X(i,j) = vc(j)*X(i,j);        
                    X(i, j) = transferFun(X(i, j), X(i, j), TFid);
                end
            end
            
            
%             % We get a new individual after the  for loop which update all
%             % decision variables of the individual.
%             % Then we bring the phase 2 of HHO to consider if there is a
%             % between position.
%             %%%%%%%%%%HHO
%             E0 = 2 * rand() - 1; %-1<E0<1
%             Escaping_Energy = E1 * (E0);  % escaping energy of rabbit
% 
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
            
        end
    end
    
%     %%%%%%%%%%OBL
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

%     % The Elite opposition based learning.
%     %%%%%%%%%%EOBL
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
%     %%%%%%%%%%
    
    convergence_curve(it)=Destination_fitness;
    it=it+1;
end
Time = toc;
end

% function X=BoundaryControl(X,lb,ub)
% [N,dim]=size(X);
% for i=1:N
%     for j=1:dim                
%         k=rand<rand; % you can change boundary-control strategy
%         if X(i,j)<lb(j)
%             if k, X(i,j)=lb(j); 
%             else X(i,j)=rand*(ub(j)-lb(j))+lb(j); 
%             end 
%         end        
%         if X(i,j)>ub(j)
%             if k, X(i,j)=ub(j);  
%             else
%                 X(i,j)=rand*(ub(j)-lb(j))+lb(j); 
%             end 
%         end
%     end
% end
% % return
% end

