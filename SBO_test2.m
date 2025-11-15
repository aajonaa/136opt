% Social based optimizer developed by Jona (upper social)
% 2024-3-18.
function [best_pos,Convergence_curve] = SBO(N, Max_FEs, lb, ub, dim, fobj)
    tic
    
    %% INITIALIZATION
    Convergence_curve = [];
    FEs = 0;

    current_X = initialization(N,dim,ub,lb);
    localElite_X = initialization(N, dim, ub, lb);
    best_pos = zeros(1,dim);
    
    if numel(lb) == 1
        lb = lb * ones(1,dim);
        ub = ub * ones(1,dim);
    end
    
    current_Fitness = inf * ones(N, 1);
    localElite_Fitness = inf * ones(N, 1);
    social_Fitness = inf * ones(N, 1);
    bestFitness = inf;
    
    for i = 1:N
        current_Fitness(i, 1) = fobj(current_X(i, :));
        FEs = FEs + 1;
        fitness = fobj(localElite_X(i ,:));
        FEs = FEs + 1;
        if current_Fitness(i, 1) < fitness
            localElite_X(i, :) = current_X(i, :);
        end
    end
    
    %% Sort the localElite fitness and best_pos
    [sorted_localElite_Fitness, idx] = sort(localElite_Fitness);
    bestFitness = sorted_localElite_Fitness(1);
    best_pos = localElite_X(idx(1), :);
    
	iter = 1;
    
    count = 0;
    
    %% Social success flag
    flag = ones(N, 1);
    
    while FEs < Max_FEs
        
        %% Select a individual from the localElite population based on the Rolette selection
        Rolette_index=RouletteWheelSelection(1./sorted_localElite_Fitness);
        if Rolette_index==-1  
            Rolette_index=1;
        end
        
        %% Update the current population
        for i = 1:N
            r1 = randn;
            r2 = randn;
            r3 = tanh(sqrt(Max_FEs - randn * FEs)/N);
%             tmp = ['complex r3 = ', num2str(r3)];
%             disp(tmp);
            r3 = real(r3);
            r4 = unifrnd(-r3, r3);
%             tmp2 = ['complex r4 = ', num2str(r4)];
%             disp(tmp2);
            if rand < r3
                for j = 1:dim
                    current_X(i, j) = (1 - r1 - r2) * current_X(i, j) + r1 * localElite_X(Rolette_index, j) + r2 * best_pos(j);
                end
            else
                count = count + 1;
                for j = 1:dim
                    current_X(i, j) = r4 * ((1 - r1 - r2) * current_X(i, j) + r1 * localElite_X(Rolette_index, j) + r2 * best_pos(j));
                end
            end
        end
        
        %% Boundary control
        current_X = BoundaryControl(current_X, lb, ub);
        
        %% Upper social: 'localElite+best不同维交换，交换位置不固定，交换数量为1' + 'localElite多维交换，交换位置固定，交换数量不固定'
        social_X = current_X;
        
        m = zeros(1, dim);
        u = randperm(dim);
        m(u(1:ceil(rand * dim))) = 1;
        for i = 1:N
            if flag(i) == 1
                social_X1 = localElite_X(i, randi(dim));
                social_X2 = best_pos(randi(dim));
                social_X(i, randi(dim)) = (social_X1 + social_X2) / 2;
            else
                for j = 1:dim
                    if m(j)
                        social_X(i, j) = localElite_X(i, j);
                    end
                end
            end
        end
        
        %% Greedy selection
        for i = 1:N
            current_Fitness(i, 1) = fobj(current_X(i, :));
            FEs = FEs + 1;
            social_Fitness(i, 1) = fobj(social_X(i, :));
            FEs = FEs + 1;
            if social_Fitness(i, 1) < current_Fitness(i, 1)
                % Social success: apply one-dimension source exchange
                flag(i, 1) = 1;
                current_X(i, :) = social_X(i, :);
            else
                % Social fail: apply multi-dimension source exchange
                flag(i, 1) = 0;
            end
        end
        
        %% Update the fitness and the localElite population
        for i = 1:N
            if current_Fitness(i, 1) < localElite_Fitness(i, 1)
                localElite_Fitness(i, 1) = current_Fitness(i, 1);
                localElite_X(i, :) = current_X(i, :);
            end
        end
        
        %% Sort the localElite fitness and best_pos
        [sorted_localElite_Fitness, idx] = sort(localElite_Fitness);
        bestFitness = sorted_localElite_Fitness(1);
        best_pos = localElite_X(idx(1), :);
    
        Convergence_curve(iter) = bestFitness;
        iter = iter + 1;
    end
    toc
    tmp3 = ['complex SBO operation2 execution count = ', num2str(count)];
    disp(tmp3);
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
 

function X = BoundaryControl(X, lb, ub)
    [N, dim] = size(X);
    for i=1:N
        for j=1:dim                
            k = rand < rand;
            if X(i,j) < lb(j)
                if k
                    X(i,j) = lb(j); 
                else
                    X(i,j) = rand*(ub(j) - lb(j)) + lb(j);
                end 
            end        
            if X(i,j) > ub(j)
                if k
                    X(i,j) = ub(j);  
                else
                    X(i,j) = rand * (ub(j) - lb(j)) + lb(j); 
                end 
            end
        end
    end
end