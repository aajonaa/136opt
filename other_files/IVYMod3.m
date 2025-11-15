function [best_pos, Convergence_curve] = IVYMod3(N, Max_FEs, lb, ub, dim, fobj)

    % Initialization
    X = zeros(N, dim);   % Initial population positions
    AllFitness = zeros(N, 1);         % Initial population costs
    GV = zeros(N, dim);         % Initial growth vectors

    FEs = 0;                    % Function evaluations counter

    for i = 1:N
        % Initialize Position (Eq. 1)
        X(i, :) = unifrnd(lb, ub, 1, dim);
        
        % Calculate Growth Vector (Eq. 6 - upper condition)
        GV(i, :) = X(i, :) ./ (ub - lb);
        
        % Evaluation
        AllFitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
    end

    % Initialize Best Cost History
    BestCosts = zeros(Max_FEs, 1);
    Convergence_curve = [];  % Preallocate for speed

    % Ivy Main Loop
    it = 1;
    while FEs <= Max_FEs
        % Initialize new population arrays
        newX = [];
        newFitness = [];
        newGV = [];
        new_pos = zeros(1, dim);

        for i = 1:N
            ii = i + 1;
            if i == N
                ii = 1;
            end

            beta_1 = 1 + (rand / 2); % beta value in Algorithm 1 (line 8)

            if AllFitness(i) < beta_1 * AllFitness(1)
                % Eq. 5-6
                new_pos = X(i, :) + abs(randn(1, dim)) .* (X(ii, :) - X(i, :)) + randn(1, dim) .* GV(i, :);
            else
                % Eq. 7
                new_pos = X(1, :) .* (rand + randn(1, dim) .* GV(i, :));
            end
            
            % Update Growth Vector (Eq. 3)
            GV(i, :) = GV(i, :) .* ((rand ^ 2) * randn(1, dim));

            % Boundary control
            new_pos = max(new_pos, lb);
            new_pos = min(new_pos, ub);

            % Update new Growth Vector (Eq. 8)
            new_gv = new_pos ./ (ub - lb);

            % Evaluate new solution
            new_cost = fobj(new_pos);
            FEs = FEs + 1;

            % Append new solutions to the population
            newX = [newX; new_pos];
            newFitness = [newFitness; new_cost];
            newGV = [newGV; new_gv];
        end

        % Merge Populations
        X = [X; newX];
        AllFitness = [AllFitness; newFitness];
        GV = [GV; newGV];

        % Sort Population by Cost
        [~, SortOrder] = sort(AllFitness);
        
        % Ensure that SortOrder size doesn't exceed array bounds
        SortOrder = SortOrder(1:N);  % Only keep the first N sorted indices

        X = X(SortOrder, :);
        AllFitness = AllFitness(SortOrder);
        GV = GV(SortOrder, :);

        % Store Best Solution Ever Found
        BestCosts(FEs) = AllFitness(1);
        Convergence_curve(it) = AllFitness(1);
        it = it + 1;
    end

    % Results
    best_pos = X(1, :);
end
