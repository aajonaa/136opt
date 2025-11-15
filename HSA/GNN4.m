function [Best_P, Convergence_curve] = GNN4(N, MaxFEs, lb, ub, dim, fobj)
%% Initialize Parameters
Best_P = zeros(1,dim);
Best_fit = inf;
AllFitness = inf*ones(N,1);
Convergence_curve = [];
it = 1;
FEs = 0;

% GNN-inspired parameters
theta = 0.6;  % Adaptive connection threshold
adj_matrix = zeros(N, N);  % Adjacency matrix
node_features = zeros(N, dim);  % Node feature matrix
learning_rate = 0.2;  % Learning rate for feature updates

%% Initialize Population
P = initialization(N, dim, ub, lb);
for i = 1:N
    FEs = FEs + 1;
    AllFitness(i) = fobj(P(i,:));
    if AllFitness(i) < Best_fit
        Best_P = P(i,:);
        Best_fit = AllFitness(i);
    end
    % Initialize node features based on position
    node_features(i,:) = (P(i,:) - lb) ./ (ub - lb);
end

%% Main Loop
while FEs <= MaxFEs
    progress_ratio = FEs / MaxFEs;
    
    %% 1. Dynamic Graph Construction
    [sorted_fitness, rank] = sort(AllFitness);
    best_id = rank(1);
    
    % Efficient similarity computation using vectorization
    distances = pdist2(P, P) ./ sqrt(dim);
    fitness_diff = pdist2(AllFitness, AllFitness) ./ (max(AllFitness) - min(AllFitness) + eps);
    
    % Compute similarities using vectorized operations
    similarity_matrix = 1 - (0.7 * distances + 0.3 * fitness_diff);
    adj_matrix = similarity_matrix >= theta;
    adj_matrix = adj_matrix - diag(diag(adj_matrix));  % Remove self-connections
    
    %% 2. Message Passing and Feature Update
    % Aggregate neighbor features
    degree = sum(adj_matrix, 2);
    neighbor_features = (adj_matrix * node_features) ./ (degree + eps);
    
    % Update node features using skip connection
    node_features = (1 - learning_rate) * node_features + learning_rate * neighbor_features;
    
    %% 3. Position Update using Graph Structure
    NewP = P;
    for i = 1:N
        neighbors = find(adj_matrix(i,:));
        if ~isempty(neighbors)
            % Local feature-based search
            local_best_idx = neighbors(argmin(AllFitness(neighbors)));
            
            % Compute attention weights
            attention = softmax(-AllFitness(neighbors));
            
            % Generate new position using attention mechanism
            if rand < 0.5  % Exploitation
                % Feature-guided local search
                direction = sum(attention .* (P(neighbors,:) - P(i,:)), 1);
                step_size = 0.1 * (1 - progress_ratio);
                NewP(i,:) = P(i,:) + step_size * direction;
            else  % Exploration
                % Graph-based global search
                global_direction = Best_P - P(i,:);
                local_direction = P(local_best_idx,:) - P(i,:);
                r1 = rand; r2 = rand;
                NewP(i,:) = P(i,:) + r1 * global_direction + r2 * local_direction;
            end
            
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
                
                % Update node features
                node_features(i,:) = (P(i,:) - lb) ./ (ub - lb);
            end
        end
    end
    
    %% 4. Adaptive Parameter Update
    if mod(it, 10) == 0
        % Adapt connection threshold
        theta = 0.8 - 0.3 * progress_ratio;
        
        % Population diversity enhancement
        if std(AllFitness) < 1e-6
            worst_indices = rank(end-floor(N*0.1)+1:end);
            for idx = worst_indices'
                P(idx,:) = lb + (ub - lb) .* rand(1, dim);
                AllFitness(idx) = fobj(P(idx,:));
                FEs = FEs + 1;
                node_features(idx,:) = (P(idx,:) - lb) ./ (ub - lb);
            end
        end
    end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

%% Helper Functions
function s = softmax(x)
    ex = exp(x - max(x));
    s = ex ./ (sum(ex) + eps);
end

function idx = argmin(x)
    [~, idx] = min(x);
end

function Positions = initialization(SearchAgents_no, dim, ub, lb)
    if length(ub) == 1
        Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
    else
        Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
    end
end