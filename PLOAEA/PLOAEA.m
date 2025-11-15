%% Polar Lights Optimization with Adaptive Elite Archive
% Polar Lights Optimization framework with the Adaptive Penalty-based Distribution strategy and the elite-driven archive mechanisms

%% APD: Adaptive PBI Distance Strategy
% Adaptive distance-based multi-objective particle swarm optimization algorithm with simple position update
% Swarm and Evolutionary Computation 94 (2025) 101890 Available online 24 February 2025

%% EDES: Elite-Driven Evolutionary Strategy
% EABC-AS: Elite-driven artificial bee colony algorithm with adaptive population scaling
% Swarm and Evolutionary Computation 94 (2025) 101893 Available online 4 March 2025
function [best_pos, Convergence_curve] = PLOAEA(N, MaxFEs, lb, ub, dim, fobj)
    %% Initialization
    FEs = 0;
    it = 1;
    AllFitness = inf * ones(N, 1);
    newFitness = inf * ones(N, 1);

    X = initialization(N, dim, ub, lb);
    V = ones(N, dim);
    newX = zeros(N, dim);

    % FIFO External Archive
    EA = 1.3; % Ratio from EABC-AS
    ArchiveSize = ceil(N * EA);
    Archive = zeros(ArchiveSize, dim); % Position storage
    ArchiveFitness = inf * ones(ArchiveSize, 1); % Fitness storage
    ArchivePtr = 1; % Pointer for FIFO

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

        % Aurora oval walk
        for i = 1:N
            a = rand() / 2 + 1;
            V(i, :) = 1 * exp((1 - a) / 100 * FEs);
            LS = V(i, :);
            GS = Levy(dim) .* (X_mean - X(i, :) + (lb + rand(1, dim) * (ub - lb)) / 2);
            newX(i, :) = X(i, :) + (w1 * LS + w2 * GS) .* rand(1, dim);
        end

        % Particle collision
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
        newX = min(max(newX, lb), ub);

        % Evaluate new solutions
        for i = 1:N
            newFitness(i) = fobj(newX(i, :));
            FEs = FEs + 1;
        end

        % Apply APD adaptively
        progress = FEs / MaxFEs;
        diversity = mean(pdist(X));
        div_threshold = 0.1 * norm(ub - lb) * (1 - progress);

        if (progress > 0.5) && (diversity < div_threshold) && (rand < 0.5)
            alpha = 0.5 - 0.3 * progress;
            beta = 0.2 + 0.8 * progress;
            subset_size = ceil(N / 4);
            worst_indices = N - subset_size + 1:N;
            subset = X(worst_indices, :);
            newX(worst_indices, :) = APD(subset, Bestpos, fobj, alpha, beta, dim);

            newX(worst_indices, :) = min(max(newX(worst_indices, :), lb), ub);
            for i = worst_indices
                newFitness(i) = fobj(newX(i, :));
                FEs = FEs + 1;
            end
        end

        % Elite-Driven Evolutionary Strategy
        Fs = 0.5 + 0.5 * (FEs / MaxFEs); % Dynamic scaling factor
        elite_size = ceil(N / 3); % Middle 33% for balance
        elite_indices = ceil(N / 3) + 1 : 2 * ceil(N / 3);
        if (progress > 0.2) && (diversity > div_threshold * 0.5) && (rand < 0.7)
            for i = elite_indices
                % Select Xpr1 from population or archive
                if rand < 0.5 && ArchivePtr > 1
                    idx = randi(min(ArchivePtr - 1, ArchiveSize));
                    Xpr1 = Archive(idx, :);
                else
                    Xpr1 = X(randi(N), :);
                end
                % Eq. 27 adaptation
                newX(i, :) = X(i, :) + Fs * (Xpr1 - X(i, :)) + Fs * (Bestpos - X(i, :));
            end
            newX(elite_indices, :) = min(max(newX(elite_indices, :), lb), ub);
            for i = elite_indices
                newFitness(i) = fobj(newX(i, :));
                FEs = FEs + 1;
            end
        end

        % Update population and archive
        for i = 1:N
            if newFitness(i) < AllFitness(i)
                % Add old solution to archive
                Archive(ArchivePtr, :) = X(i, :);
                ArchiveFitness(ArchivePtr) = AllFitness(i);
                ArchivePtr = mod(ArchivePtr, ArchiveSize) + 1;
                if ArchivePtr == 1
                    ArchivePtr = ArchiveSize;
                end
                % Update population
                X(i, :) = newX(i, :);
                AllFitness(i) = newFitness(i);
            end
        end

        % Sort and update best solution
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

% Revised APD function
function updated_population = APD(population, best_solution, objective_func, alpha, beta, dim)
    [num_solutions, ~] = size(population);
    updated_population = zeros(num_solutions, dim);
    f_best = objective_func(best_solution);

    for i = 1:num_solutions
        candidate = population(i, :);
        f_x = objective_func(candidate);
        lambda_vec = randn(1, dim);
        lambda_vec = lambda_vec / norm(lambda_vec);

        d1 = dot((candidate - best_solution), lambda_vec);
        d2 = norm(candidate - (best_solution + d1 * lambda_vec));

        improvement = max(0, (f_x - f_best) / (abs(f_best) + eps));
        update_scale = min(1, improvement);
        update_direction = -alpha * d1 * lambda_vec + beta * d2 * randn(1, dim);
        new_candidate = candidate + update_scale * update_direction;

        updated_population(i, :) = new_candidate;
    end
end

% Levy and initialization functions
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1 / beta);
    o = step;
end