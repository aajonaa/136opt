% binary version 4.15
function [best_pos,Convergence_curve]=RMRIME(N,Max_FEs,lb,ub,dim,fobj)
    % disp('RIME is now tackling your problem')

    
    lb = 0;
    ub = 1;
    
    % initialize position
    best_pos=zeros(1,dim);
    bestFitness=inf;%change this to -inf for maximization problems
%     X=initialization(N,dim,ub,lb);%Initialize the set of random solutions
    X=initialization(N,dim,ub,lb) > 0.5;%Initialize the set of random solutions
    X2 = X;
    Lb=lb.*ones(1,dim);% lower boundary 
    Ub=ub.*ones(1,dim);% upper boundary
    % it=1;%Number of iterations
%     FEs=0;
    iter = 1;
    Convergence_curve=[];
    AllFitness=zeros(1,N);%Initialize the fitness value
    newAllFitness=zeros(1,N);
    W = 5;%Soft-rime parameters, discussed in subsection 4.3.1 of the paper
    %Calculate the fitness value of the initial position
    for i=1:N
        AllFitness(1,i)=fobj(X(i,:));%Calculate the fitness value for each search agent
%         FEs=FEs+1;
        %Make greedy selections
        if AllFitness(1,i)<bestFitness
            bestFitness=AllFitness(1,i);
            best_pos=X(i,:);
        end
    end
    % Main loop
    while iter < Max_FEs

        [sorted_AllFitness, ~] = sort(AllFitness);
        Rolette_index=RouletteWheelSelection(1./sorted_AllFitness);
        if Rolette_index==-1  
            Rolette_index=1;
        end

        F = (rand-0.5)*2*cos((pi*iter/(Max_FEs/10)))*(1-round(iter*W/Max_FEs)/W);%Parameters of Eq.(3),(4),(5)
        p = sqrt(iter/Max_FEs);%Eq.(6)
        newX = X;%Recording new populations
        normalizedAllFitness=normr(AllFitness);%Parameters of Eq.(7)
        for i=1:N
            for j=1:dim
                %Soft-rime search strategy
                r1=rand();
                if r1< p
    %                 newX(i,j)=best_pos(1,j)+ F *((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
%                     newX(i,j) = X(Rolette_index, j) + F *((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
                    temp = X(Rolette_index, j) + F *((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
                    newX(i,j) = transferFun(newX(i, j), temp, TFid);
                end
                %Hard-rime puncture mechanism
                r2=rand();
                if r2 < normalizedAllFitness(i)
    %                 newX(i,j)=best_pos(1,j);%Eq.(7)
%                     newX(i,j) = X(Rolette_index, j);%Eq.(7)
                    temp = X(Rolette_index, j);%Eq.(7)
                    newX(i,j) = transferFun(newX(i, j), temp, TFid);
                end
            end
        end
        for i=1:N
            %Boundary absorption
%             Flag4ub=newX(i,:)>ub;
%             Flag4lb=newX(i,:)<lb;
%             newX(i,:)=(newX(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            newAllFitness(1,i)=fobj(newX(i,:));
%             FEs=FEs+1;
            %Positive greedy selection mechanism
            if newAllFitness(1,i)<AllFitness(1,i)
                AllFitness(1,i) = newAllFitness(1,i);
                X(i,:) = newX(i,:);
                if newAllFitness(1,i)< bestFitness
                   bestFitness=AllFitness(1,i);
                   best_pos=X(i,:);
                end
            end
        end
        
        % crossover
        history_X = X;
        for i = 1:N
            m = zeros(1, dim);
            u = randperm(dim);
            m(1, u(1:floor(rand*dim+1))) = 1;

    %         X2(i, :) = m .* X(i, :) + randn * ~m .* (X(floor(rand * N + 1), :)- X(i, :));
%             X2(i, :) = history_X(i, :) + randn * m .* (history_X(Rolette_index, :) - history_X(i, :));
            temp = history_X(i, :) + randn * m .* (history_X(Rolette_index, :) - history_X(i, :));
            X2(i, :) = transferFun(X2(i, :), temp, TFid);
        end

        for i=1:N

            % Return back the search agents that go beyond the boundaries of the search space
%             Flag4ub=X2(i,:)>ub;
%             Flag4lb=X2(i,:)<lb;
%             X2(i,:)=(X2(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

            % Calculate objective function for each search agent
            fitness = fobj(X2(i, :));
%             FEs = FEs + 1;

            if fitness < AllFitness(1, i) 
                AllFitness(1, i) = fitness;
                X(i, :) = X2(i, :);
            end

        %             if fitness<bestFitness % Change this to > for maximization problem
            if AllFitness(1, i)<bestFitness % Change this to > for maximization problem
        %                 bestFitness=fitness; % Update alpha
                bestFitness=AllFitness(1, i); % Update alpha
                best_pos=X(i,:);
            end
        end
        
        Convergence_curve(iter)=bestFitness;
        iter=iter+1;
    end
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


