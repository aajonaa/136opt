function [Best_P, Convergence_curve] = GNN2(N, MaxFEs, lb, ub, dim, fobj)
%% Initialize Parameters
Best_P = zeros(1,dim);
Best_fit = inf;
AllFitness = inf*ones(N,1);
Convergence_curve = [];
FEs = 0;
it = 1;

% GNN-inspired parameters
hidden_dim = 32;  % Hidden dimension for node features
num_layers = 3;   % Number of message passing layers
learning_rate = 0.1; % Learning rate for feature updates

% Initialize node features (analogous to GNN node embeddings)
node_features = randn(N, hidden_dim);  % Random initial node features
edge_threshold = 0.5; % Threshold for edge creation

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
    % Current progress ratio
    progress_ratio = FEs / MaxFEs;
    
    %% 1. Graph Construction Phase (Message Passing Analogue)
    % Construct adjacency matrix using node features
    adj_matrix = construct_graph(node_features, edge_threshold);
    
    % Update node features through message passing
    for layer = 1:num_layers
        node_features = message_passing(node_features, adj_matrix, P, AllFitness);
    end
    
    %% 2. Node Feature Aggregation and Update
    for i = 1:N
        % Aggregate neighboring features
        neighbor_idx = find(adj_matrix(i,:));
        if ~isempty(neighbor_idx)
            aggregated_features = mean(node_features(neighbor_idx,:), 1);
            
            % Update node features using aggregated information
            node_features(i,:) = node_features(i,:) + ...
                learning_rate * (aggregated_features - node_features(i,:));
            
            % Generate new position based on aggregated features
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
    
    %% 3. Topology Adaptation (Graph Structure Learning)
    if mod(it, 10) == 0  % Adapt topology periodically
        edge_threshold = adapt_threshold(progress_ratio, edge_threshold);
        node_features = update_features(node_features, P, AllFitness, Best_P);
    end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

% Helper Functions
function adj_matrix = construct_graph(node_features, threshold)
    % Construct adjacency matrix using cosine similarity of node features
    N = size(node_features, 1);
    similarity_matrix = zeros(N, N);
    
    for i = 1:N
        for j = i+1:N
            similarity = cosine_similarity(node_features(i,:), node_features(j,:));
            if similarity > threshold
                similarity_matrix(i,j) = similarity;
                similarity_matrix(j,i) = similarity;
            end
        end
    end
    
    adj_matrix = similarity_matrix > threshold;
end

function new_features = message_passing(features, adj_matrix, positions, fitness)
    % Implement GNN-style message passing
    N = size(features, 1);
    new_features = features;
    
    for i = 1:N
        neighbors = find(adj_matrix(i,:));
        if ~isempty(neighbors)
            % Aggregate neighbor features weighted by fitness
            neighbor_weights = softmax(-fitness(neighbors));
            weighted_features = features(neighbors,:) .* neighbor_weights';
            new_features(i,:) = mean(weighted_features, 1);
        end
    end
end

function new_position = generate_position(node_feature, current_pos, best_pos, lb, ub)
    % Generate new position using node features as guidance
    feature_direction = tanh(node_feature(1:length(current_pos)));
    step_size = 0.1 * (1 + rand);
    
    % Combine current position, best position, and feature direction
    new_position = current_pos + ...
        step_size * feature_direction .* (best_pos - current_pos);
    
    % Bound constraints
    new_position = min(max(new_position, lb), ub);
end

function new_threshold = adapt_threshold(progress_ratio, current_threshold)
    % Adapt edge threshold based on optimization progress
    min_threshold = 0.3;
    max_threshold = 0.7;
    new_threshold = max_threshold - progress_ratio * (max_threshold - min_threshold);
end

function new_features = update_features(features, positions, fitness, best_pos)
    % Update node features based on optimization state
    N = size(features, 1);
    dim = size(features, 2);
    
    % Calculate relative performance
    normalized_fitness = (fitness - min(fitness)) / (max(fitness) - min(fitness) + eps);
    
    % Update features based on optimization performance
    for i = 1:N
        performance_factor = exp(-normalized_fitness(i));
        position_difference = positions(i,:) - best_pos;
        
        % Combine current features with position information
        new_features(i,:) = features(i,:) + ...
            0.1 * performance_factor * randn(1, dim);
    end
end

function similarity = cosine_similarity(v1, v2)
    % Compute cosine similarity between two vectors
    similarity = dot(v1, v2) / (norm(v1) * norm(v2) + eps);
end

function weights = softmax(x)
    % Compute softmax of input vector
    exp_x = exp(x - max(x));
    weights = exp_x / sum(exp_x);
end