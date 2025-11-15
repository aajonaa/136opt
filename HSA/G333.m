% %% Part 1: Graph Initialization and Structure
% function [best_pos, Convergence_curve] = G333(N, MaxFEs, lb, ub, dim, fobj)
%     % Initialize basic parameters
%     best_pos = zeros(1, dim);
%     bestFitness = inf;
%     AllFitness = inf * ones(N, 1);
%     Convergence_curve = [];
%     FEs = 0;
%     it = 1;
% 
%     % Enhanced GNN parameters
%     n_heads = 3;                  % Number of attention heads
%     embedding_dim = dim;          % Feature embedding dimension
%     hidden_dim = round(dim * 1.5); % Hidden layer dimension
% 
%     % Initialize population and node features
%     P = initialization(N, dim, lb, ub);
%     node_features = zeros(N, embedding_dim);
% 
%     % Initial fitness evaluation
%     for i = 1:N
%         FEs = FEs + 1;
%         AllFitness(i) = fobj(P(i,:));
%         % Initialize node features with position and fitness information
%         node_features(i,:) = create_node_features(P(i,:), AllFitness(i), lb, ub);
%         if AllFitness(i) < bestFitness
%             best_pos = P(i,:);
%             bestFitness = AllFitness(i);
%         end
%     end
% 
%     %% Part 2: GNN Message Passing and Update
%     while FEs <= MaxFEs
%         progress_ratio = FEs / MaxFEs;
% 
%         % Multi-head attention mechanism
%         attention_weights = zeros(N, N, n_heads);
%         for head = 1:n_heads
%             attention_weights(:,:,head) = compute_attention(node_features, progress_ratio);
%         end
% 
%         % Aggregate messages using attention
%         messages = aggregate_messages(node_features, attention_weights, hidden_dim);
% 
%         % Update node features through neural network layers
%         updated_features = update_node_features(messages, node_features, hidden_dim);
% 
%         % Convert updated features back to solution space
%         new_positions = decode_features_to_positions(updated_features, lb, ub);
% 
%         %% Part 3: Solution Update and Diversity Management
%         for i = 1:N
%             % Evaluate new position
%             new_fitness = fobj(new_positions(i,:));
%             FEs = FEs + 1;
% 
%             % Update if improved
%             if new_fitness < AllFitness(i)
%                 P(i,:) = new_positions(i,:);
%                 AllFitness(i) = new_fitness;
%                 node_features(i,:) = create_node_features(P(i,:), new_fitness, lb, ub);
% 
%                 if new_fitness < bestFitness
%                     best_pos = P(i,:);
%                     bestFitness = new_fitness;
%                 end
%             end
%         end
% 
%         % Diversity enhancement through feature space
%         if mod(it, 10) == 0
%             P = enhance_diversity(P, node_features, AllFitness, lb, ub);
%         end
% 
%         it = it + 1;
%         Convergence_curve(it) = bestFitness;
%     end
% end
% 
% %% Helper Functions
% function features = create_node_features(position, fitness, lb, ub)
%     % Create richer node features combining position and fitness information
%     normalized_pos = (position - lb) ./ (ub - lb);
%     features = [normalized_pos, fitness];
% end
% 
% function attention = compute_attention(features, progress_ratio)
%     % Compute attention scores using scaled dot-product attention
%     scaling_factor = sqrt(size(features, 2));
%     scores = (features * features') / scaling_factor;
% 
%     % Apply temperature scaling based on progress
%     temperature = 1 + progress_ratio;
%     scores = scores / temperature;
% 
%     % Softmax normalization
%     attention = exp(scores) ./ sum(exp(scores), 2);
% end
% 
% function messages = aggregate_messages(features, attention_weights, hidden_dim)
%     % Aggregate messages using multi-head attention
%     n_heads = size(attention_weights, 3);
%     messages = zeros(size(features, 1), hidden_dim);
% 
%     for head = 1:n_heads
%         head_message = attention_weights(:,:,head) * features;
%         messages = messages + head_message;
%     end
%     messages = messages / n_heads;
% end
% 
% function updated_features = update_node_features(messages, features, hidden_dim)
%     % Update features using a simple feed-forward network
%     W1 = randn(size(features, 2), hidden_dim) / sqrt(hidden_dim);
%     W2 = randn(hidden_dim, size(features, 2)) / sqrt(size(features, 2));
% 
%     % Forward pass with ReLU activation
%     hidden = max(0, messages * W1);
%     updated_features = hidden * W2;
% 
%     % Residual connection
%     updated_features = updated_features + features;
% end
% 
% function positions = decode_features_to_positions(features, lb, ub)
%     % Convert features back to position space
%     positions = features(:, 1:length(lb)) .* (ub - lb) + lb;
% end
% 
% function P_new = enhance_diversity(P, features, fitness, lb, ub)
%     % Enhance diversity through feature space perturbation
%     [~, worst_idx] = maxk(fitness, round(size(P, 1) * 0.1));
% 
%     for i = worst_idx'
%         % Generate new solution in feature space
%         P_new = P;
%         P_new(i,:) = lb + (ub - lb) .* rand(1, size(P, 2));
%     end
% end

function [best_pos, Convergence_curve] = EnhancedGNN(N, MaxFEs, lb, ub, dim, fobj)
    % Initialize basic parameters
    best_pos = zeros(1, dim);
    bestFitness = inf;
    AllFitness = inf * ones(N, 1);
    Convergence_curve = [];
    FEs = 0;
    it = 1;

    % Enhanced GNN parameters
    n_heads = 3;                  % Number of attention heads
    embedding_dim = dim;          % Feature embedding dimension - now matches problem dimension
    hidden_dim = round(dim * 1.5); % Hidden layer dimension
    
    % Initialize population and node features
    P = initialization(N, dim, lb, ub);
    node_features = zeros(N, embedding_dim);  % Now matches problem dimension
    
    % Initial fitness evaluation
    for i = 1:N
        FEs = FEs + 1;
        AllFitness(i) = fobj(P(i,:));
        % Initialize node features - now only using normalized position
        node_features(i,:) = create_node_features(P(i,:), lb, ub);
        if AllFitness(i) < bestFitness
            best_pos = P(i,:);
            bestFitness = AllFitness(i);
        end
    end

    while FEs <= MaxFEs
        progress_ratio = FEs / MaxFEs;
        
        % Multi-head attention mechanism
        attention_weights = zeros(N, N, n_heads);
        for head = 1:n_heads
            % Now passing AllFitness separately for attention computation
            attention_weights(:,:,head) = compute_attention(node_features, AllFitness, progress_ratio);
        end
        
        % Aggregate messages using attention
        messages = aggregate_messages(node_features, attention_weights, hidden_dim);
        
        % Update node features through neural network layers
        updated_features = update_node_features(messages, node_features, hidden_dim);
        
        % Convert updated features back to solution space
        new_positions = decode_features_to_positions(updated_features, lb, ub);
        
        % Solution Update and Diversity Management
        for i = 1:N
            new_fitness = fobj(new_positions(i,:));
            FEs = FEs + 1;
            
            if new_fitness < AllFitness(i)
                P(i,:) = new_positions(i,:);
                AllFitness(i) = new_fitness;
                node_features(i,:) = create_node_features(P(i,:), lb, ub);
                
                if new_fitness < bestFitness
                    best_pos = P(i,:);
                    bestFitness = new_fitness;
                end
            end
        end
        
        % Diversity enhancement
        if mod(it, 10) == 0
            P = enhance_diversity(P, node_features, AllFitness, lb, ub);
        end
        
        it = it + 1;
        Convergence_curve(it) = bestFitness;
    end
end

%% Helper Functions
function features = create_node_features(position, lb, ub)
    % Create node features using only normalized position
    features = (position - lb) ./ (ub - lb);
end

function attention = compute_attention(features, fitness, progress_ratio)
    % Compute attention scores using both position features and fitness
    N = size(features, 1);
    
    % Compute position-based similarity
    scaling_factor = sqrt(size(features, 2));
    position_scores = (features * features') / scaling_factor;
    
    % Compute fitness-based similarity
    fitness_diff = abs(fitness - fitness') / (max(fitness) - min(fitness));
    fitness_scores = 1 ./ (1 + fitness_diff);
    
    % Combine position and fitness scores
    scores = (1 - progress_ratio) * position_scores + progress_ratio * fitness_scores;
    
    % Apply temperature scaling
    temperature = 1 + progress_ratio;
    scores = scores / temperature;
    
    % Softmax normalization
    attention = exp(scores) ./ sum(exp(scores), 2);
end

function messages = aggregate_messages(features, attention_weights, hidden_dim)
    n_heads = size(attention_weights, 3);
    messages = zeros(size(features, 1), hidden_dim);
    
    for head = 1:n_heads
        % Project features to hidden dimension first
        W = randn(size(features, 2), hidden_dim) / sqrt(hidden_dim);
        head_message = attention_weights(:,:,head) * (features * W);
        messages = messages + head_message;
    end
    messages = messages / n_heads;
end

function updated_features = update_node_features(messages, features, hidden_dim)
    % Project original features to hidden dimension
    W1 = randn(size(features, 2), hidden_dim) / sqrt(hidden_dim);
    W2 = randn(hidden_dim, size(features, 2)) / sqrt(size(features, 2));
    
    % Forward pass with ReLU activation
    hidden = max(0, messages + (features * W1));  % Added skip connection
    updated_features = hidden * W2;
    
    % Residual connection
    updated_features = updated_features + features;
end

function positions = decode_features_to_positions(features, lb, ub)
    % Convert features back to position space
    positions = features .* (ub - lb) + lb;
end

function P_new = enhance_diversity(P, features, fitness, lb, ub)
    P_new = P;
    [~, worst_idx] = maxk(fitness, round(size(P, 1) * 0.1));
    
    for i = worst_idx'
        P_new(i,:) = lb + (ub - lb) .* rand(1, size(P, 2));
    end
end