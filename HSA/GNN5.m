function [Best_P, Convergence_curve] = GNN5(N, MaxFEs, lb, ub, dim, fobj)
%% Initialize Parameters
Best_P = zeros(1, dim);
Best_fit = inf;
AllFitness = inf * ones(N, 1);
Convergence_curve = [];
FEs = 0;
it = 1;

% Simple parameters
connection_threshold = 0.5;  % Fixed threshold for connections
learning_rate = 0.3;        % Fixed learning rate

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
    %% 1. Simple Graph Construction
    % Build connections based on fitness and distance
    connections = zeros(N, N);
    for i = 1:N
        for j = i+1:N
            % Simple similarity based on normalized distance and fitness difference
            dist = norm(P(i,:) - P(j,:)) / norm(ub - lb);
            fit_diff = abs(AllFitness(i) - AllFitness(j)) / (max(AllFitness) - min(AllFitness));
            
            % Combined similarity measure
            similarity = 1 - (0.5 * dist + 0.5 * fit_diff);
            
            % Establish connection if similarity is high enough
            if similarity > connection_threshold
                connections(i,j) = 1;
                connections(j,i) = 1;
            end
        end
    end
    
    %% 2. Position Update
    NewP = P;
    for i = 1:N
        % Find connected neighbors
        neighbors = find(connections(i,:));
        
        if ~isempty(neighbors)
            % Simple weighted average of neighbor positions
            neighbor_weights = 1 ./ (AllFitness(neighbors) + eps);
            weights = neighbor_weights / sum(neighbor_weights);
            
            % Generate new position
            mean_pos = sum(P(neighbors,:) .* weights', 1);
            
            % Move towards mean position and best solution
            r1 = rand; r2 = rand;
            NewP(i,:) = P(i,:) + ...
                       learning_rate * (r1 * (mean_pos - P(i,:)) + ...
                                     r2 * (Best_P - P(i,:)));
            
            % Boundary handling
            NewP(i,:) = min(max(NewP(i,:), lb), ub);
            
            % Evaluate new position
            new_fitness = fobj(NewP(i,:));
            FEs = FEs + 1;
            
            % Update if better
            if new_fitness < AllFitness(i)
                P(i,:) = NewP(i,:);
                AllFitness(i) = new_fitness;
                if new_fitness < Best_fit
                    Best_P = NewP(i,:);
                    Best_fit = new_fitness;
                end
            end
        end
    end
    
    %% 3. Simple Diversity Maintenance
    if mod(it, 20) == 0  % Every 20 iterations
        % Reinitialize worst solution
        [~, worst_idx] = max(AllFitness);
        P(worst_idx,:) = lb + (ub - lb) .* rand(1, dim);
        AllFitness(worst_idx) = fobj(P(worst_idx,:));
        FEs = FEs + 1;
    end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

%% Basic initialization function
function Positions = initialization(SearchAgents_no, dim, ub, lb)
    if length(ub) == 1
        Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
    else
        Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
    end
end