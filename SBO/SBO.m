%% Status-based Optimization: Algorithm and comprehensive performance analysis (SBO)
function [best_pos,Convergence_curve] = SBO(N, Max_FEs, lb, ub, dim, fobj)
    
    %INITIALIZATION
    Convergence_curve = [];
    FEs = 0;


    current_X = initialization(N,dim,ub,lb); % see Eq.1 in [1]
    localElite_X = initialization(N, dim, ub, lb);
        
    current_Fitness = inf * ones(N, 1);
    localElite_Fitness = inf * ones(N, 1);
    social_Fitness = inf * ones(N, 1);
    
    for i = 1:N
        current_Fitness(i, 1) = fobj(current_X(i, :));
        FEs = FEs + 1;
        fitness = fobj(localElite_X(i ,:));
        FEs = FEs + 1;
        if current_Fitness(i, 1) < fitness
            localElite_X(i, :) = current_X(i, :);
            localElite_Fitness(i, 1) = current_Fitness(i, 1);
        else
            localElite_Fitness(i, 1) = fitness;
        end
    end
    
    %% Sort the best_pos
    [sorted_localElite_Fitness, idx] = sort(localElite_Fitness);
    best_pos = localElite_X(idx(1), :);
    bestFitness = sorted_localElite_Fitness(1);
    
	iter = 1;
    
    %% Social success flag
    flag = ones(N, 1);
    
    while FEs < Max_FEs
        
        %% Select a individual from the localElite population based on the Roulette selection
        Roulette_index=RouletteWheelSelection(1./(sorted_localElite_Fitness + eps));
        if Roulette_index==-1  
            Roulette_index=1;
        end
        
        %% Update the current population
        for i = 1:N
            w1 = randn;
            w2 = randn;
            w3 = tanh((sqrt(abs(Max_FEs - randn * FEs))/i)^(FEs/Max_FEs));
            w4 = unifrnd(-w3, w3);
            if rand < w3
                for j = 1:dim
                    current_X(i, j) = (1 - w1 - w2) * current_X(i, j) + w1 * localElite_X(Roulette_index, j) + w2 * best_pos(j);
                end
            else
                for j = 1:dim
                    current_X(i, j) = w4 * ((1 - w1 - w2) * current_X(i, j) + w1 * localElite_X(Roulette_index, j) + w2 * best_pos(j));
                end
            end
        end
        
        %% Boundary control
        current_X = BoundaryControl(current_X, lb, ub);
        
        %% Upward social
        social_X = current_X;
        
            for i = 1:N
                if flag(i) == 1
                    social_X1 = localElite_X(i, randi(dim));
                    social_X2 = best_pos(randi(dim));
                    social_X(i, randi(dim)) = (social_X1 + social_X2) / 2;
                end
            end
            
            m = zeros(1, dim);
            u = randperm(dim);
            m(u(1:ceil(rand * dim))) = 1;
            for i = 1:N
                if flag(i) == 0
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
                current_Fitness(i, 1) = social_Fitness(i, 1);
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
        if sorted_localElite_Fitness(1) < bestFitness
            bestFitness = sorted_localElite_Fitness(1);
            best_pos = localElite_X(idx(1), :);
        end
    
        Convergence_curve(iter) = bestFitness;
        iter = iter + 1;
    end
end  

function X = BoundaryControl(X, low, up)
    [N, dim] = size(X);
    
    if isscalar(low)
        low = repmat(low, 1, dim);
    end
    if isscalar(up)
        up = repmat(up, 1, dim);
    end
    
    for i = 1:N
        for j = 1:dim                
            k = rand < rand; 
            
            if X(i,j) < low(j) 
                if k
                    X(i,j) = low(j); 
                else
                    X(i,j) = rand * (up(j) - low(j)) + low(j);  
                end 
            end        
            
            if X(i,j) > up(j)  
                if k
                    X(i,j) = up(j);  
                else
                    X(i,j) = rand * (up(j) - low(j)) + low(j); 
                end 
            end
        end
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

function X = initialization(N, dim, ub, lb)
    if isscalar(lb)
        lb = repmat(lb, 1, dim);
    end
    if isscalar(ub)
        ub = repmat(ub, 1, dim);
    end
    X = rand(N, dim) .* (ub - lb) + lb;
end