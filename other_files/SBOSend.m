% Social based optimizer developed by Jona (upward social)
% 2024-3-18.
function [best_pos,Convergence_curve] = SBO(N, Max_FEs, lb, ub, dim, fobj)
    tic
    
    %% INITIALIZATION
    Convergence_curve = [];
    FEs = 0;

    X = initialization(N,dim,ub,lb);
    XE = initialization(N, dim, ub, lb);
    best_pos = zeros(1,dim);
    
    if numel(lb) == 1
        lb = lb * ones(1,dim);
        ub = ub * ones(1,dim);
    end
    
    Fitness = inf * ones(N, 1);
    FitnessE = inf * ones(N, 1);
    FitnessS = inf * ones(N, 1);
    bestFitness = inf;
    
    for i = 1:N
        Fitness(i, 1) = fobj(X(i, :));
        FEs = FEs + 1;
        fitness = fobj(XE(i ,:));
        FEs = FEs + 1;
        if Fitness(i, 1) < fitness
            XE(i, :) = X(i, :);
        end
    end
    
    %% Sort the localElite fitness and best_pos
    [sorted_FitnessE, idx] = sort(FitnessE);
    best_pos = XE(idx(1), :);
    
	iter = 1;

    %% Social success flag
    flag = ones(N, 1);
    
    while FEs < Max_FEs
        
        %% Select a individual from the localElite population based on the Rolette Wheel selection
        Rolette_index=RouletteWheelSelection(1./sorted_FitnessE);
        if Rolette_index==-1  
            Rolette_index=1;
        end
        
        %% Update the current population
        for i = 1:N
            R1 = randn;
            R2 = randn;
            R3 = tanh((sqrt(abs(Max_FEs - randn * FEs))/i)^(FEs/Max_FEs));
            R4 = real(unifrnd(-R3, R3));
            if rand < R3
                for j = 1:dim
                    X(i, j) = (1 - R1 - R2) * X(i, j) + R1 * XE(Rolette_index, j) + R2 * best_pos(j);
                end
            else
                for j = 1:dim
                    X(i, j) = R4 * ((1 - R1 - R2) * X(i, j) + R1 * XE(Rolette_index, j) + R2 * best_pos(j));
                end
            end
        end
        
        %% Boundary control
        X = BoundaryControl(X, lb, ub);
        
        %% Upword social: 'localElite+best不同维交换，交换位置不固定，交换数量为1' + 'localElite多维交换，交换位置固定，交换数量不固定'
        XS = X;
        
        m = zeros(1, dim);
        u = randperm(dim);
        m(u(1:ceil(rand * dim))) = 1;
        for i = 1:N
            if flag(i) == 1
                XS1 = XE(i, randi(dim));
                XS2 = best_pos(randi(dim));
                XS(i, randi(dim)) = (XS1 + XS2) / 2;
            else
                for j = 1:dim
                    if m(j)
                        XS(i, j) = XE(i, j);
                    end
                end
            end
        end
        
        %% Greedy selection
        for i = 1:N
            Fitness(i, 1) = fobj(X(i, :));
            FEs = FEs + 1;
            FitnessS(i, 1) = fobj(XS(i, :));
            FEs = FEs + 1;
            if FitnessS(i, 1) < Fitness(i, 1)
                % Social success: apply one-dimension source exchange
                flag(i, 1) = 1;
                X(i, :) = XS(i, :);
            else
                % Social fail: apply multi-dimension source exchange
                flag(i, 1) = 0;
            end
        end
        
        %% Update the fitness and the localElite population
        for i = 1:N
            if Fitness(i, 1) < FitnessE(i, 1)
                FitnessE(i, 1) = Fitness(i, 1);
                XE(i, :) = X(i, :);
            end
        end
        
        %% Sort the localElite fitness and best_pos
        [sorted_FitnessE, idx] = sort(FitnessE);
        bestFitness = sorted_FitnessE(1);
        best_pos = XE(idx(1), :);
    
        Convergence_curve(iter) = bestFitness;
        iter = iter + 1;
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