% function [best_pos, Convergence_curve] = APDPLO(N, MaxFEs, lb, ub, dim, fobj)
%     %% Initialization
%     FEs = 0;
%     it = 1;
%     AllFitness = inf * ones(N, 1);
%     newFitness = inf * ones(N, 1);
% 
%     X = initialization(N, dim, ub, lb);
%     V = ones(N, dim);
%     newX = zeros(N, dim);
% 
%     for i = 1:N
%         AllFitness(i) = fobj(X(i, :));
%         FEs = FEs + 1;
%     end
% 
%     [AllFitness, SortOrder] = sort(AllFitness);
%     X = X(SortOrder, :);
%     Bestpos = X(1, :);
%     bestFitness = AllFitness(1);
% 
%     Convergence_curve = [];
%     Convergence_curve(it) = bestFitness;
% 
%     %% Main loop
%     while FEs <= MaxFEs
%         X_sum = sum(X, 1);
%         X_mean = X_sum / N;
%         w1 = tansig((FEs / MaxFEs)^4);
%         w2 = exp(-(2 * FEs / MaxFEs)^3);
% 
%         % Aurora oval walk (PLO core)
%         for i = 1:N
%             a = rand() / 2 + 1;
%             V(i, :) = 1 * exp((1 - a) / 100 * FEs);
%             LS = V(i, :);
%             GS = Levy(dim) .* (X_mean - X(i, :) + (lb + rand(1, dim) * (ub - lb)) / 2);
%             newX(i, :) = X(i, :) + (w1 * LS + w2 * GS) .* rand(1, dim);
%         end
% 
%         % Particle collision (PLO core)
%         E = sqrt(FEs / MaxFEs);
%         A = randperm(N);
%         for i = 1:N
%             for j = 1:dim
%                 if (rand < 0.05) && (rand < E)
%                     newX(i, j) = X(i, j) + sin(rand * pi) * (X(i, j) - X(A(i), j));
%                 end
%             end
%         end
% 
%         % Boundary control
%         Flag4ub = newX > ub;
%         Flag4lb = newX < lb;
%         newX = (newX .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
% 
%         % Evaluate new solutions
%         for i = 1:N
%             newFitness(i) = fobj(newX(i, :));
%             FEs = FEs + 1;
%         end
% 
%         % Apply APD adaptively to the worst 25% -- REVISED
%         progress = FEs / MaxFEs;
%         diversity = mean(pdist(X));
%         div_threshold = 0.1 * norm(ub - lb);
% 
%         if (progress > 0.5) && (diversity < div_threshold) && (rand < 0.5) % Changed to progress > 0.5
%             alpha = 0.5 - 0.3 * progress; % Decreases from 0.5 to 0.2
%             beta = 0.2 + 0.8 * progress; % Increases from 0.2 to 1.0
%             subset_size = ceil(N / 4);
%             worst_indices = N - subset_size + 1 : N; % Worst 25%
%             subset = X(worst_indices, :);
%             newX(worst_indices, :) = APD(subset, Bestpos, fobj, alpha, beta, dim); % Updated APD call
% 
%             % Reapply boundary control
%             Flag4ub = newX(worst_indices, :) > ub;
%             Flag4lb = newX(worst_indices, :) < lb;
%             newX(worst_indices, :) = (newX(worst_indices, :) .* (~(Flag4ub + Flag4lb))) + ...
%                                      ub .* Flag4ub + lb .* Flag4lb;
% 
%             % Re-evaluate after APD
%             for i = worst_indices
%                 newFitness(i) = fobj(newX(i, :));
%                 FEs = FEs + 1;
%             end
%         end
% 
%         % Update population
%         for i = 1:N
%             if newFitness(i) < AllFitness(i)
%                 X(i, :) = newX(i, :);
%                 AllFitness(i) = newFitness(i);
%             end
%         end
% 
%         % Sort and update the best solution
%         [AllFitness, SortOrder] = sort(AllFitness);
%         X = X(SortOrder, :);
%         if AllFitness(1) < bestFitness
%             Bestpos = X(1, :);
%             bestFitness = AllFitness(1);
%         end
% 
%         it = it + 1;
%         Convergence_curve(it) = bestFitness;
%         best_pos = Bestpos;
%     end
% end
% 
% % Revised APD function
% function updated_population = APD(population, best_solution, objective_func, alpha, beta, dim)
%     [num_solutions, ~] = size(population);
%     updated_population = zeros(num_solutions, dim);
%     f_best = objective_func(best_solution);
% 
%     for i = 1:num_solutions
%         candidate = population(i, :);
%         f_x = objective_func(candidate);
%         % Unique lambda_vec for each solution -- REVISED
%         lambda_vec = randn(1, dim);
%         lambda_vec = lambda_vec / norm(lambda_vec);
% 
%         d1 = dot((candidate - best_solution), lambda_vec);
%         d2 = norm(candidate - (best_solution + d1 * lambda_vec));
% 
%         improvement = max(0, (f_x - f_best) / (abs(f_best) + eps));
%         update_scale = min(1, improvement);
%         update_direction = -alpha * d1 * lambda_vec + beta * d2 * randn(1, dim);
%         new_candidate = candidate + update_scale * update_direction;
% 
%         updated_population(i, :) = new_candidate;
%     end
% end
% 
% % Levy and initialization functions remain unchanged
% function o = Levy(d)
%     beta = 1.5;
%     sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
%     u = randn(1, d) * sigma;
%     v = randn(1, d);
%     step = u ./ abs(v).^(1 / beta);
%     o = step;
% end

function [best_pos, Convergence_curve] = APDPLO(N, MaxFEs, lb, ub, dim, fobj)
    %% Initialization
    FEs = 0;
    it = 1;
    AllFitness = inf * ones(N, 1);
    newFitness = inf * ones(N, 1);

    X = initialization(N, dim, ub, lb);
    V = ones(N, dim);
    newX = zeros(N, dim);

    for i = 1:N
        AllFitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
    end

    [AllFitness, SortOrder] = sort(AllFitness);
    X = X(SortOrder, :);
    Bestpos = X(1, :);
    bestFitness = AllFitness(1);

    Convergence_curve = [];
    Convergence_curve(it) = bestFitness;

    %% Main loop
    while FEs <= MaxFEs
        X_sum = sum(X, 1);
        X_mean = X_sum / N;
        w1 = tansig((FEs / MaxFEs)^4);
        w2 = exp(-(2 * FEs / MaxFEs)^3);

        % Aurora oval walk (PLO core)
        for i = 1:N
            a = rand() / 2 + 1;
            V(i, :) = 1 * exp((1 - a) / 100 * FEs);
            LS = V(i, :);
            GS = Levy(dim) .* (X_mean - X(i, :) + (lb + rand(1, dim) * (ub - lb)) / 2);
            newX(i, :) = X(i, :) + (w1 * LS + w2 * GS) .* rand(1, dim);
        end

        % Particle collision (PLO core)
        E = sqrt(FEs / MaxFEs);
        A = randperm(N);
        for i = 1:N
            for j = 1:dim
                if (rand < 0.05) && (rand < E)
                    newX(i, j) = X(i, j) + sin(rand * pi) * (X(i, j) - X(A(i), j));
                end
            end
        end

        % Boundary control
        Flag4ub = newX > ub;
        Flag4lb = newX < lb;
        newX = (newX .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

        % Evaluate new solutions
        for i = 1:N
            newFitness(i) = fobj(newX(i, :));
            FEs = FEs + 1;
        end

        % Apply APD adaptively with improved strategy
        progress = FEs / MaxFEs;
        diversity = mean(pdist(X)); 
        scale_factor = 0.05 * norm(ub - lb); 
        div_threshold = max(scale_factor, 0.2 * norm(ub - lb) * (1 - progress)); 

        if (progress > 0.2) && (diversity < div_threshold) && (rand < 0.8) 
            alpha = 0.15 + 0.25 * (1 - progress); 
            beta = 0.3 + 0.4 * progress; 
            lambda_vec = ones(1, dim) + 0.03 * randn(1, dim); 
            lambda_vec = lambda_vec / norm(lambda_vec);

            subset_size = ceil(N / 3); 
            selected_idx = randperm(N, subset_size);
            newX(selected_idx, :) = APD(newX(selected_idx, :), Bestpos, fobj, lambda_vec, alpha, beta);

            % Reapply boundary control after APD
            Flag4ub = newX(selected_idx, :) > ub;
            Flag4lb = newX(selected_idx, :) < lb;
            newX(selected_idx, :) = (newX(selected_idx, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            % Re-evaluate after APD
            for i = 1:subset_size
                newFitness(selected_idx(i)) = fobj(newX(selected_idx(i), :));
                FEs = FEs + 1;
            end
        end

        % Update population
        for i = 1:N
            if newFitness(i) < AllFitness(i)
                X(i, :) = newX(i, :);
                AllFitness(i) = newFitness(i);
            end
        end

        % Sort and update the best solution
        [AllFitness, SortOrder] = sort(AllFitness);
        X = X(SortOrder, :);
        if AllFitness(1) < bestFitness
            Bestpos = X(1, :);
            bestFitness = AllFitness(1);
        end

        it = it + 1;
        Convergence_curve(it) = bestFitness;
        best_pos = Bestpos;
    end
end

function updated_population = APD(population, best_solution, objective_func, lambda_vec, alpha, beta)
    [num_solutions, num_variables] = size(population);
    updated_population = zeros(num_solutions, num_variables);
    f_best = objective_func(best_solution);

    for i = 1:num_solutions
        candidate = population(i, :);
        f_x = objective_func(candidate);

        % Compute d1 and d2 for improved PBI
        d1 = dot((candidate - best_solution), lambda_vec);
        d2 = norm(candidate - (best_solution + d1 * lambda_vec));

        % Enhanced adaptive update with fitness-distance correlation
        fitness_ratio = (f_x - f_best) / (abs(f_best) + eps);
        improvement = max(0, fitness_ratio);
        distance_factor = exp(-d2 / (1e-6 + norm(ub - lb)));
        
        update_scale = min(1, 0.5 + 0.5 * improvement) * distance_factor;
        direction = -alpha * d1 * lambda_vec + beta * d2 * (randn(1, num_variables) * 0.5 + 0.5);
        new_candidate = candidate + update_scale * direction;

        % Dynamic bounds for exploration
        new_candidate = new_candidate + 0.1 * Levy(num_variables) .* (ub - lb) .* (1 - update_scale);

        updated_population(i, :) = new_candidate;
    end
end

% Levy and initialization functions remain unchanged