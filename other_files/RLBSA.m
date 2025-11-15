%% Claude 2025-2-15
function [best_pos, Convergence_curve] = BSA(N, Max_FEs, lb, ub, dim, fobj)
    tic
    Convergence_curve = [];
    FEs = 0;
    bestFitness = inf;
    best_pos = zeros(1, dim);

    %% Initialize Python environment if needed
    if count(py.sys.path, '') == 0
        insert(py.sys.path, int32(0), '');
    end

    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    end

    AllFitness = inf * ones(1, N);
    X = initialization(N, dim, ub, lb);
    historyX = initialization(N, dim, ub, lb);
    it = 1;

    % Calculate the fitness of the initial population
    for i = 1:N
        AllFitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
        if AllFitness(i) < bestFitness
            bestFitness = AllFitness(i);
            best_pos = X(i, :);
        end
    end

    % Main loop: Selection, Mutation, and Crossover
    while FEs < Max_FEs
        if rand < rand
            historyX = X;
        end

        historyX = historyX(randperm(N), :);

        %% 🔹 **Use RL Attention for Selection**
        try
            % Convert MATLAB array to Python numpy array
            py_population = py.numpy.array(X);

            % Call the Python function
            py_result = py.rl_attention.get_attention_modified_population(py_population);

            % Convert Python result back to MATLAB array
            newX_cell = cell(py_result);
            newX = cellfun(@double, newX_cell, 'UniformOutput', false);
            newX = cell2mat(newX);

            % Reshape if necessary to match original dimensions
            if size(newX, 1) ~= N || size(newX, 2) ~= dim
                newX = reshape(newX, N, dim);
            end
        catch ME
            warning('Python interface failed: %s. Using default mutation.', ME.message);
            % Fallback strategy: simple mutation
            newX = X + 0.1 * randn(size(X));
        end

        %% Boundaries Control
        newX = BoundaryControl(newX, lb, ub);

        %% Selection-II
        newAllFitness = zeros(1, N);  % Initialize newAllFitness
        for i = 1:N
            newAllFitness(i) = fobj(newX(i, :));
            FEs = FEs + 1;
        end

        ind = newAllFitness < AllFitness;
        AllFitness(ind) = newAllFitness(ind);
        X(ind, :) = newX(ind, :);

        [globalminimum, ind] = min(AllFitness);
        best_pos = X(ind, :);
        bestFitness = globalminimum;
        Convergence_curve(it) = bestFitness;
        it = it + 1;
           
        % Optional: Print progress every 10 iterations
        if mod(it, 10) == 0
            fprintf('Iteration: %d, Best Fitness: %.4e\n', it, bestFitness);
        end 

    end
    toc
end

function X = BoundaryControl(X, lb, ub)
    % Apply boundary constraints
    X = max(min(X, ub), lb);
end


% F=get_scale_factor; % see Eq.5 in [1], you can define other F generation strategies 
% map=zeros(N,dim); % see Algorithm-2 in [1]  
% if rand<rand
%     for i=1:N  
%         u=randperm(dim); 
%         map(i,u(1:ceil(DIM_RATE*rand*dim)))=1;
%     end
% else
%     for i=1:N  
%         map(i,randi(dim))=1;
%     end
% end
% newX=X+(map.*F).*(historyX-X);   % see Eq.5 in [1]  