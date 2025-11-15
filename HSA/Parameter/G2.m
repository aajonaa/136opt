function [Best_P, Convergence_curve] = GNN_Fast3(N, MaxFEs, lb, ub, dim, fobj)
%% This algorithm implements a Graph Neural Network (GNN) inspired metaheuristic
% It combines concepts from GNNs with optimization techniques:
% 1. Graph Structure: Solutions form nodes in a dynamic graph
% 2. Message Passing: Information flows between connected solutions
% 3. Node Features: Solutions evolve based on neighborhood information
% 4. Attention Mechanism: Adaptive edge weights based on solution similarity

%% Initialize Parameters
Best_P = zeros(1, dim);        % Best position found so far
Best_fit = inf;                % Best fitness value
AllFitness = inf * ones(N, 1); % Fitness values for all solutions
Convergence_curve = [];        % Track convergence history
FEs = 0;                       % Function evaluation counter
it = 1;                        % Iteration counter

% GNN-inspired parameters
% message_passing_strength = 0.3;  % Controls information flow strength (analogous to neural message passing)
% attention_weight = 0.5;         % Balances position vs fitness similarity in attention
% graph_sparsity = 4;            % Controls graph connectivity (like k-nearest neighbors in GNNs)
update_interval = 3;           % Frequency of graph structure updates
theta = 0.5;                   % Edge formation threshold (similar to attention threshold in GNNs)

% Initialize graph structure (representing the solution space as a graph)
connections = sparse(N, N);     % Adjacency matrix (sparse for efficiency)
similarity = zeros(N, N);       % Pairwise solution similarities (edge weights)
NoUpdateCount = zeros(1, N);    % Track solution stagnation
last_update_iter = 0;          % Last graph update iteration

%% Initialize Population (Graph Nodes)
P = initialization(N, dim, lb, ub);  % Each solution becomes a node in our graph
for i = 1:N
    FEs = FEs + 1;
    AllFitness(i) = fobj(P(i,:));
    if AllFitness(i) < Best_fit
        Best_P = P(i,:);
        Best_fit = AllFitness(i);
    end
end

while FEs <= MaxFEs
    progressRatio = FEs / MaxFEs;  % Used for adaptive parameter control
    
    %% 1. Graph Structure Update with Attention
    % This section implements the graph attention mechanism similar to GAT (Graph Attention Networks)
    if it - last_update_iter >= update_interval
        connections = sparse(N, N);
        [SortFit, FitRank] = sort(AllFitness);
        
        for i = 1:N
            for j = i+1:N
                % Calculate normalized distance (position-based similarity)
                NormDistIJ = sqrt(sum(((P(i,:) - P(j,:)) ./ (ub - lb)) .^ 2)) / sqrt(dim);
                
                % Calculate fitness difference (fitness-based similarity)
                DiffFitIJ = abs(AllFitness(i) - AllFitness(j)) / (max(AllFitness) - min(AllFitness));
                
                % Attention mechanism: Combine position and fitness similarities
                % Uses a logarithmic attention function inspired by transformer architectures
                similarity(i,j) = 1 - log(1 + ((1-progressRatio) * NormDistIJ + progressRatio * DiffFitIJ) * (exp(1) - 1));
                similarity(j,i) = similarity(i,j);
                
                % Edge formation with threshold (like attention dropout in GNNs)
                if similarity(i,j) >= theta
                    connections(i,j) = 1;
                    connections(j,i) = 1;
                end
            end
        end
        last_update_iter = it;
    end
    
    %% 2. Message Passing and Position Update
    % Implements neural message passing where solutions exchange information
    AAF = 0.5 * 2^(1 - progressRatio);  % Adaptive adjustment factor (like learning rate in GNNs)
    
    for i = 1:N
        neighbors = find(connections(i,:));  % Get neighborhood (receptive field)
        if ~isempty(neighbors)
            % Aggregate neighborhood information (analogous to neighborhood aggregation in GNNs)
            LocalMean = mean(P(neighbors,:), 1);
            
            % Levy flight for exploration (adds stochasticity to message passing)
            R = Levy(1,dim) .* (rand(1,dim) > 0.5);
            w1 = 1 - progressRatio^2;  % Dynamic weighting factor
            
            % Generate new solutions using three different update rules
            % These are analogous to different types of neural message passing
            
            % Rule 1: Local-global balance
            Node1 = P(i,:) + AAF * R .* (w1 * (Best_P - LocalMean) + (1-w1) * (LocalMean - P(i,:)));
            
            % Rule 2: Global influence
            Node2 = Best_P + AAF * R .* (w1 * (LocalMean - P(i,:)) + (1-w1) * (P(i,:) - Best_P));
            
            % Rule 3: Local influence
            Node3 = LocalMean + AAF * R .* (w1 * (P(i,:) - Best_P) + (1-w1) * (Best_P - LocalMean));
            
            % Boundary handling
            Node1 = min(max(Node1, lb), ub);
            Node2 = min(max(Node2, lb), ub);
            Node3 = min(max(Node3, lb), ub);
            
            % Evaluate candidate solutions
            F_Node1 = fobj(Node1);
            F_Node2 = fobj(Node2);
            F_Node3 = fobj(Node3);
            FEs = FEs + 3;
            
            % Select best candidate (like max-pooling in GNNs)
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
    
    %% 3. Population Diversity Management
    % Similar to dropout in GNNs - helps prevent overfitting to local optima
    [~, pruneIdx] = maxk(NoUpdateCount, 1 + floor(rand * N * 0.1));
    for IDX = pruneIdx
        % Reinitialize stagnant solutions (analogous to neuron reactivation)
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
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

% Levy flight implementation for enhanced exploration
% Generates steps following a Levy distribution (heavy-tailed)
function o = Levy(n,d)
    beta = 1.5;  % Levy index (controls the distribution's heavy tail)
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(n,d) * sigma;
    v = randn(n,d);
    step = u./abs(v).^(1/beta);
    o = step;
end

% Initialize population with uniform random distribution
function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Positions = lb + (ub - lb) .* rand(SearchAgents_no, dim);
end