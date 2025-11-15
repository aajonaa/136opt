function [Best_P, Convergence_curve] = GNN_Fast3(N, MaxFEs, lb, ub, dim, fobj)
%% GNN-Fast3: A Graph Neural Network Inspired Metaheuristic Algorithm
% Input Parameters:
%   N       - Population size
%   MaxFEs  - Maximum function evaluations
%   lb      - Lower bounds vector
%   ub      - Upper bounds vector
%   dim     - Problem dimension
%   fobj    - Objective function handle
% Outputs:
%   Best_P  - Best solution found
%   Convergence_curve - Convergence history

%% Initialize Parameters
Best_P = zeros(1, dim);
Best_fit = inf;
AllFitness = inf * ones(N, 1);
Convergence_curve = [];
FEs = 0;
it = 1;

% GNN core parameters with adaptive capabilities
message_passing_strength = 0.8;  % Controls information flow between solutions
attention_weight = 0.3;         % Balances position vs. fitness influence
graph_sparsity = min(N-1, max(3, floor(N/4)));  % Adaptive neighborhood size
update_interval = 3;            % Graph update frequency
theta = 0.5;                   % Connection threshold

% Initialize tracking variables
connections = sparse(N, N);     % Graph structure (adjacency matrix)
NoUpdateCount = zeros(1, N);    % Stagnation counter
last_update_iter = 0;          % Last graph update iteration
stagnation_threshold = 5;      % Threshold for diversity injection
diversity_history = zeros(1, 10); % Track population diversity
improvement_rate = 1;          % Track optimization progress

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

prev_best_fit = Best_fit;  % For improvement rate calculation

while FEs <= MaxFEs
    progressRatio = FEs / MaxFEs;
    
    % Update GNN parameters based on search progress
    message_passing_strength = adaptMessagePassingStrength(progressRatio, improvement_rate);
    attention_weight = adaptAttentionWeight(progressRatio, diversity_history(end));
    graph_sparsity = adaptGraphSparsity(N, progressRatio, improvement_rate);
    
    %% 1. Graph Structure Update with Dynamic Attention
    if it - last_update_iter >= update_interval
        connections = sparse(N, N);
        
        % Update diversity measure
        diversity = calculateDiversity(P, lb, ub);
        diversity_history = [diversity_history(2:end), diversity];
        
        % Build graph connections with attention mechanism
        for i = 1:N
            similarities_i = zeros(1, N);
            for j = (i+1):N
                pos_similarity = calculatePositionSimilarity(P(i,:), P(j,:), lb, ub);
                fit_similarity = calculateFitnessSimilarity(AllFitness(i), AllFitness(j), AllFitness);
                
                % Combine similarities using attention weight
                similarity_ij = combineNodeSimilarities(pos_similarity, fit_similarity, attention_weight);
                similarities_i(j) = similarity_ij;
                
                % Establish symmetric connections above threshold
                if similarity_ij >= theta
                    connections(i,j) = similarity_ij;
                    connections(j,i) = similarity_ij;
                end
            end
        end
        last_update_iter = it;
    end
    
    %% 2. Message Passing and Solution Update
    AAF = 0.5 * 2^(1 - progressRatio) * (1 + diversity_history(end));
    
    for i = 1:N
        neighbors = find(connections(i,:));
        if ~isempty(neighbors)
            % Compute weighted neighborhood information
            weights = connections(i,neighbors);
            weights = weights / sum(weights);
            LocalMean = calculateWeightedLocalMean(P, neighbors, weights, message_passing_strength, i);
            
            % Generate candidate solutions using Lévy flight
            R = Levy(1,dim) .* (rand(1,dim) > progressRatio);
            w1 = calculateDynamicWeight(progressRatio, improvement_rate);
            
            % Create three candidate solutions using different update rules
            [Node1, Node2, Node3] = generateCandidateSolutions(P(i,:), Best_P, LocalMean, AAF, R, w1, message_passing_strength);
            
            % Apply boundary constraints
            Node1 = min(max(Node1, lb), ub);
            Node2 = min(max(Node2, lb), ub);
            Node3 = min(max(Node3, lb), ub);
            
            % Evaluate candidates
            F_Node1 = fobj(Node1);
            F_Node2 = fobj(Node2);
            F_Node3 = fobj(Node3);
            FEs = FEs + 3;
            
            % Select best candidate
            [MinFit, idx] = min([F_Node1, F_Node2, F_Node3, AllFitness(i)]);
            if MinFit < AllFitness(i)
                switch idx
                    case 1
                        P(i,:) = Node1;
                    case 2
                        P(i,:) = Node2;
                    case 3
                        P(i,:) = Node3;
                end
                AllFitness(i) = MinFit;
                NoUpdateCount(i) = 0;
                
                if MinFit < Best_fit
                    Best_P = P(i,:);
                    Best_fit = MinFit;
                end
            else
                NoUpdateCount(i) = NoUpdateCount(i) + 1;
            end
        end
    end
    
    % Update improvement rate for adaptive parameters
    improvement_rate = (prev_best_fit - Best_fit) / (prev_best_fit + eps);
    prev_best_fit = Best_fit;
    
    %% 3. Diversity Management
    if mean(NoUpdateCount) > stagnation_threshold
        [~, pruneIdx] = maxk(NoUpdateCount, 1 + floor(0.1 * N));
        for IDX = pruneIdx
            % Reinitialize stagnant solutions
            P(IDX,:) = lb + (ub - lb) .* rand(1, dim);
            connections(IDX,:) = 0;
            connections(:,IDX) = 0;
            AllFitness(IDX) = fobj(P(IDX,:));
            FEs = FEs + 1;
            NoUpdateCount(IDX) = 0;
            
            if AllFitness(IDX) < Best_fit
                Best_P = P(IDX,:);
                Best_fit = AllFitness(IDX);
            end
        end
    end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

%% Helper Functions for GNN Operations
function strength = adaptMessagePassingStrength(progress, improvement)
    % Adapts message passing strength based on search progress and improvement rate
    base_strength = 0.8 * (1 - progress^2);
    strength = base_strength * (1 + improvement);
    strength = min(max(strength, 0.2), 0.8);
end

function weight = adaptAttentionWeight(progress, diversity)
    % Adapts attention weight based on search progress and population diversity
    weight = 0.3 + 0.4 * progress * (1 - diversity);
    weight = min(max(weight, 0.2), 0.8);
end

function sparsity = adaptGraphSparsity(N, progress, improvement)
    % Adapts graph connectivity based on population size and search state
    base_sparsity = floor(N/4);
    if improvement < 0.001
        base_sparsity = base_sparsity + 2;
    end
    sparsity = min(N-1, max(3, floor(base_sparsity * (1 - 0.5 * progress))));
end

function diversity = calculateDiversity(P, lb, ub)
    % Calculates normalized population diversity
    ranges = ub - lb;
    std_pos = std(P);
    diversity = mean(std_pos ./ ranges);
end

function sim = calculatePositionSimilarity(pos1, pos2, lb, ub)
    % Calculates normalized position-based similarity between solutions
    ranges = ub - lb;
    norm_diff = (pos1 - pos2) ./ ranges;
    sim = 1 / (1 + sqrt(mean(norm_diff.^2)));
end

function sim = calculateFitnessSimilarity(fit1, fit2, all_fits)
    % Calculates normalized fitness-based similarity between solutions
    range = max(all_fits) - min(all_fits) + eps;
    sim = 1 - abs(fit1 - fit2) / range;
end

function combined = combineNodeSimilarities(pos_sim, fit_sim, attention_weight)
    % Combines position and fitness similarities using attention mechanism
    combined = attention_weight * pos_sim + (1 - attention_weight) * fit_sim;
end

function mean_pos = calculateWeightedLocalMean(P, neighbors, weights, message_strength, current)
    % Calculates weighted mean of neighborhood with message passing strength
    weighted_sum = weights * P(neighbors,:);
    mean_pos = message_strength * weighted_sum + (1 - message_strength) * P(current,:);
end

function w = calculateDynamicWeight(progress, improvement)
    % Calculates dynamic weight factor based on search progress
    w = (1 - progress^2) * (1 + improvement);
    w = min(max(w, 0.1), 0.9);
end

function [n1, n2, n3] = generateCandidateSolutions(current, best, local, AAF, R, w1, msg_strength)
    % Generates three candidate solutions using different update strategies
    n1 = current + AAF * R .* (w1 * (best - local) + (1-w1) * msg_strength * (local - current));
    n2 = best + AAF * R .* (w1 * msg_strength * (local - current) + (1-w1) * (current - best));
    n3 = local + AAF * R .* (w1 * (current - best) + (1-w1) * msg_strength * (best - local));
end

function o = Levy(n,d)
    % Generates Lévy flight step sizes for enhanced exploration
    beta = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(n,d) * sigma;
    v = randn(n,d);
    o = u./abs(v).^(1/beta);
end

function Positions = initialization(SearchAgents_no, dim, ub, lb)
    % Initializes population uniformly within bounds
    Positions = lb + (ub - lb) .* rand(SearchAgents_no, dim);
end