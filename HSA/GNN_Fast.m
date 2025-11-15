function [Best_P, Convergence_curve] = GNN_low_performance(N, MaxFEs, lb, ub, dim, fobj)
%% Initialize Parameters
Best_P = zeros(1, dim);
Best_fit = inf;
AllFitness = inf * ones(N, 1);
Convergence_curve = [];
FEs = 0;
it = 1;

% Efficiency parameters
update_prob = 0.3;     % Probability of updating each solution (like dropout)
k_neighbors = 3;       % Maximum number of neighbors to consider (sparsity)
update_interval = 5;   % Update graph every n iterations

% Initialize sparse connection matrix
connections = sparse(N, N);
last_update_iter = 0;

%% Initialize Population
P = initialization(N, dim, lb, ub);
for i = 1:N
    FEs = FEs + 1;
    AllFitness(i) = fobj(P(i,:));
    if AllFitness(i) < Best_fit
        Best_P = P(i,:);
        Best_fit = AllFitness(i);
    end
end

while FEs <= MaxFEs
    %% 1. Selective Graph Update (every update_interval iterations)
    if it - last_update_iter >= update_interval
        % Sort population by fitness
        [~, sorted_idx] = sort(AllFitness);
        
        % Update connections using k-nearest neighbors approach
        connections = sparse(N, N);
        for i = 1:N
            % Only connect to k best solutions among nearby solutions
            start_idx = max(1, i - k_neighbors);
            end_idx = min(N, i + k_neighbors);
            potential_neighbors = sorted_idx(start_idx:end_idx);
            
            % Add connections
            for j = potential_neighbors'
                if i ~= j
                    connections(i,j) = 1;
                    connections(j,i) = 1;
                end
            end
        end
        
        last_update_iter = it;
    end
    
    %% 2. Fast Position Update with Dropout
    % Randomly select solutions to update
    update_mask = rand(N, 1) < update_prob;
    update_indices = find(update_mask);
    
    % Update selected solutions in parallel
    for i = update_indices'
        neighbors = find(connections(i,:));
        
        if ~isempty(neighbors)
            % Quick neighbor selection (take only the best ones)
            [~, best_neighbor_idx] = min(AllFitness(neighbors));
            best_neighbor = neighbors(best_neighbor_idx);
            
            % Simple position update using only the best neighbor
            r = rand;
            NewP = P(i,:) + r * (P(best_neighbor,:) - P(i,:)) + ...
                   (1-r) * (Best_P - P(i,:));
            
            % Quick boundary handling
            NewP = min(max(NewP, lb), ub);
            
            % Evaluate new position
            new_fitness = fobj(NewP);
            FEs = FEs + 1;
            
            % Update if better
            if new_fitness < AllFitness(i)
                P(i,:) = NewP;
                AllFitness(i) = new_fitness;
                if new_fitness < Best_fit
                    Best_P = NewP;
                    Best_fit = new_fitness;
                end
            end
        end
    end
    
    %% 3. Quick Diversity Check
    if mod(it, 20) == 0 && FEs < MaxFEs * 0.8  % Only in early and middle stages
        % Find worst solution
        [~, worst_idx] = max(AllFitness);
        
        % Simple reset with influence from best solution
        P(worst_idx,:) = 0.5 * (lb + (ub - lb) .* rand(1, dim)) + ...
                         0.5 * Best_P;
        AllFitness(worst_idx) = fobj(P(worst_idx,:));
        FEs = FEs + 1;
    end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

%% Simple initialization
function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Positions = lb + (ub - lb) .* rand(SearchAgents_no, dim);
end