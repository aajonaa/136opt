function [Best_P, Convergence_curve] = GNN_MA(N, MaxFEs, lb, ub, dim, fobj)
%% Initialize Parameters
Best_P = zeros(1,dim);
Best_fit = inf;
AllFitness = inf*ones(N,1);
Convergence_curve = [];
FEs = 0;
it = 1;

% GNN-inspired parameters (reduced complexity)
hidden_dim = dim;  % Match hidden dimension with problem dimension for simplicity
edge_threshold = 0.5; % Threshold for edge creation

% Initialize node features matching the problem dimension
node_features = randn(N, hidden_dim);

%% Initialize Population (Nodes)
P = initialization(N, dim, ub, lb);
for i = 1:N
    FEs = FEs + 1;
    AllFitness(i) = fobj(P(i,:));
    if AllFitness(i) < Best_fit
        Best_P = P(i,:);
        Best_fit = AllFitness(i);
    end
end

while FEs <= MaxFEs
    progress_ratio = FEs / MaxFEs;
    
    %% 1. Graph Construction Phase
    adj_matrix = construct_graph(node_features, edge_threshold);
    
    %% 2. Node Feature Update and Position Generation
    for i = 1:N
        % Get neighbor indices
        neighbors = find(adj_matrix(i,:));
        
        if ~isempty(neighbors)
            % Update node features based on neighbors
            neighbor_features = node_features(neighbors,:);
            neighbor_fitness = AllFitness(neighbors);
            
            % Calculate weights based on fitness
            weights = softmax(-neighbor_fitness);
            weights = weights(:);  % Ensure column vector
            
            % Weighted aggregation of neighbor features
            aggregated_features = zeros(1, hidden_dim);
            for j = 1:length(neighbors)
                aggregated_features = aggregated_features + weights(j) * neighbor_features(j,:);
            end
            
            % Update current node features
            node_features(i,:) = 0.7 * node_features(i,:) + 0.3 * aggregated_features;
            
            % Generate new position
            new_position = generate_position(node_features(i,:), P(i,:), Best_P, lb, ub);
            
            % Evaluate new position
            new_fitness = fobj(new_position);
            FEs = FEs + 1;
            
            % Update if better
            if new_fitness < AllFitness(i)
                P(i,:) = new_position;
                AllFitness(i) = new_fitness;
                if new_fitness < Best_fit
                    Best_P = new_position;
                    Best_fit = new_fitness;
                end
            end
        end
    end
    
    %% 3. Topology Adaptation
    if mod(it, 10) == 0
        edge_threshold = adapt_threshold(progress_ratio);
        node_features = update_features(node_features, P, AllFitness, Best_P);
    end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

function adj_matrix = construct_graph(node_features, threshold)
    N = size(node_features, 1);
    adj_matrix = zeros(N, N);
    
    for i = 1:N
        for j = i+1:N
            % Compute cosine similarity
            similarity = max(min(dot(node_features(i,:), node_features(j,:)) / ...
                (norm(node_features(i,:)) * norm(node_features(j,:)) + eps), 1), -1);
            
            if similarity > threshold
                adj_matrix(i,j) = 1;
                adj_matrix(j,i) = 1;
            end
        end
    end
end

function new_position = generate_position(node_feature, current_pos, best_pos, lb, ub)
    % Use node features to guide the search direction
    direction = tanh(node_feature);  % Squash values to [-1, 1]
    
    % Random step size
    step_size = 0.1 * (1 + rand);
    
    % Generate new position
    new_position = current_pos + step_size * direction .* (best_pos - current_pos);
    
    % Bound constraints
    new_position = min(max(new_position, lb), ub);
end

function new_threshold = adapt_threshold(progress_ratio)
    % Simple linear adaptation
    min_threshold = 0.3;
    max_threshold = 0.7;
    new_threshold = max_threshold - progress_ratio * (max_threshold - min_threshold);
end

function new_features = update_features(features, positions, fitness, best_pos)
    N = size(features, 1);
    dim = size(features, 2);
    new_features = features;
    
    % Normalize fitness to [0,1]
    normalized_fitness = (fitness - min(fitness)) / (max(fitness) - min(fitness) + eps);
    
    for i = 1:N
        % Calculate influence factor based on fitness
        influence = exp(-normalized_fitness(i));
        
        % Update features with small random perturbation
        new_features(i,:) = features(i,:) + 0.1 * influence * randn(1, dim);
    end
end

function weights = softmax(x)
    % Numerically stable softmax
    x = x - max(x);  % For numerical stability
    exp_x = exp(x);
    weights = exp_x / (sum(exp_x) + eps);
end

% Basic initialization function if not provided
function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Boundary_no = size(ub, 2);
    if Boundary_no == 1
        Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
    else
        for i = 1:dim
            ub_i = ub(i);
            lb_i = lb(i);
            Positions(:, i) = rand(SearchAgents_no, 1) .* (ub_i - lb_i) + lb_i;
        end
    end
end