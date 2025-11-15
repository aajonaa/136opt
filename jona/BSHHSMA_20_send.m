% Slime Mold Algorithm modified by Jona 2023-10-29.
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
function [best_pos, convergence_curve] = BSHHSMA_20_send(N, Max_FEs, lb, ub, dim, fobj)
tic
disp('BSHHSMA_20_send is now tackling your problem')

% Initialize position
best_pos = zeros(1, dim);
Destination_fitness = inf; % Change this to -inf for maximization problems
AllFitness = inf * ones(N, 1); % Record the fitness of all slime mold
weight = ones(N, dim); % Fitness weight of each slime mold
% Initialize the set of random solutions
X = initialization(N, dim, ub, lb); 
convergence_curve = [];
it = 1;  % Number of iterations
lb = ones(1, dim) .* lb; % Lower boundary 
ub = ones(1, dim) .* ub; % Upper boundary
z = 0.03; % Parameter

%%%%%%%%%% BSA
% Initialize the initial BSA population and historical population.
for i = 1:N
    bs_X(i, :) = X(i, :);
end
bs_h_X = X;
%%%%%%%%%%

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

    [sorted_fitness_value, sorted_fitness_index] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = sorted_fitness_value(N);
    bestFitness = sorted_fitness_value(1);
    
    %%%%%%%%%% BSA
    % Insert the BSA to SMA before 1/3 iteration for extensive exploration.
    if FEs < Max_FEs / 3
        % Initialize the BS(backtracking search) population's fitness to inf.
        bs_X_fitness = inf * ones(1, N);
        % Use for loop to assigned the SMA's population and fitness
        % value to the counterpart of the BS(backtracking search).(fitness
        % increasing).
        j = 1;
        for i = 1:N
            bs_X_fitness(j) = sorted_fitness_value(i, :);
            bs_X(j, :) = X(sorted_fitness_index(i), :);
            j = j + 1;
        end
        % Randomly get the historical backtracking search population by
        % permutating the rows of the backtracking search population.
        bs_h_X = bs_X(randperm(N), :);
        % Scale facotr used in map.
        F = 3 * randn;
        % Initialize the map(logical matrix) to be 0.
        map = zeros(N, dim);

        % Randomly set the decision variables to be 1 used for mutaion.
        for i = 1:30
            u = randperm(dim);
            map(i, u(1:ceil(rand * dim))) = 1;
        end
        
        % Get the next generation by mutation.
        offsprings = bs_X + (map .* F) .* (bs_h_X - bs_X);
        % Control the boundaries of the decision variables.
        offsprings = BoundaryControl(offsprings, lb, ub);

        % Initialize the offsprings fitness to be inf.
        fitness_offsprings = inf * ones(1, N);

        % Caculate the offsprings fitness.
        for i = 1:N
            fitness_offsprings(i) = fobj(offsprings(i, :));
            FEs = FEs + 1;
        end

        % To choose the best N individuals from two generation.
        bs_X_fitness = reshape(bs_X_fitness, [1, N]); % I add this in order to avoid error when running the program.
        
        % To get the logical array for select the best fitness.
        ind = fitness_offsprings < bs_X_fitness;
        
        % Combined the X and offsprings to form a better population X.
        bs_X_fitness(ind) = fitness_offsprings(ind);
        bs_X(ind, :) = offsprings(ind, :);

        % We need to assign the population's fitness value to the SMA's
        % AllFitness after update the population.
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
    %%%%%%%%%% end BSA

    S = bestFitness - worstFitness + eps;  % Plus eps to avoid denominator zero

    % Calculate the fitness weight of each slime mold
    for i = 1:N
        for j = 1:dim
            if i <= (N / 2)  %Eq.(2.5)
                weight(sorted_fitness_index(i), j) = 1 + rand() * log10((bestFitness - sorted_fitness_value(i)) / S + 1);
            else
                weight(sorted_fitness_index(i), j) = 1 - rand() * log10((bestFitness - sorted_fitness_value(i)) / S + 1);
            end
        end
    end
    
    % Update the best fitness value and best position
    if bestFitness < Destination_fitness
        best_pos = X(sorted_fitness_index(1), :);
        Destination_fitness = bestFitness;
    end
    
    a = atanh( -(FEs / Max_FEs) + 1);   %Eq.(2.4)
    b = 1 - FEs / Max_FEs;
    
    % Update the Position of search agents
    for i = 1:N
        % Random distribution status.
        if rand < z     %Eq.(2.7)
            X(i, :) = (ub - lb) * rand + lb;
        % Exploration status.
        else
            p = tanh(abs(AllFitness(i) - Destination_fitness));  %Eq.(2.2) P is between [0, 1]. Decreasing from 1 to 0. 1-->0
            vb = unifrnd(-a, a, 1, dim);  %Eq.(2.3)
            vc = unifrnd(-b, b, 1, dim);
            for j = 1:dim
                r = rand();
                AA = randi([1, N]);  % Two positions randomly selected from population
                B = randi([1, N]);
                % Exploration status.
                if r < p    %Eq.(2.1)
                    X(i, j) = best_pos(j) + vb(j) * (weight(i, j) * X(AA,j) - X(B, j));
                    
                % Exploitation status.
                else
                    X(i, j) = vc(j) * X(i, j);                    
                end
            end 
        end
    end    
    convergence_curve(it) = Destination_fitness;
    it = it + 1;
end
toc
end

% BSA boundary control.
function X = BoundaryControl(X, lb, ub)
[N, dim] = size(X);
for i = 1:N
    for j = 1:dim                
        k = rand < rand; 
        if X(i, j)<lb(j)
            if k, X(i, j)=lb(j); 
            else X(i, j) = rand * (ub(j) - lb(j)) + lb(j); 
            end 
        end        
        if X(i, j) > ub(j)
            if k, X(i, j) = ub(j);  
            else
                X(i, j) = rand * (ub(j) - lb(j)) + lb(j); 
            end 
        end
    end
end
end

