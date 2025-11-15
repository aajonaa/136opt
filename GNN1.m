function [Best_P, Convergence_curve] = GNN1(N, MaxFEs, lb, ub, dim, fobj)
%% Initialize Parameters
Best_P = zeros(1, dim);
Best_fit = inf;
AllFitness = inf * ones(N, 1);
Convergence_curve = [];
FEs = 0;
it = 1;

% Enhanced efficiency parameters inspired by GNN concepts
message_passing_strength = 0.3;  % Controls information flow between nodes
attention_weight = 0.5;         % Attention mechanism weight
graph_sparsity = 4;            % Number of neighbors (like sparse attention)
update_interval = 3;           % Graph update frequency

% Initialize adaptive parameters (inspired by MGO)
beta = 0.8;                    % Controlled exploration factor
gamma = 1/sqrt(1-power(beta,2)); % Relativistic factor
w = 0.8;                       % Inertia weight

% Initialize graph structure
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
    %% 1. Dynamic Graph Structure Update
    if it - last_update_iter >= update_interval
        % Calculate normalized fitness for attention mechanism
        norm_fitness = (AllFitness - min(AllFitness)) / (max(AllFitness) - min(AllFitness));
        
        % Update graph using attention-based connectivity
        connections = sparse(N, N);
        [~, sorted_idx] = sort(AllFitness);
        
        for i = 1:N
            % Attention-based neighbor selection
            attention_scores = exp(-norm_fitness/attention_weight);
            [~, neighbor_idx] = maxk(attention_scores, graph_sparsity);
            
            % Establish weighted connections
            for j = neighbor_idx'
                if i ~= j
                    % Edge weight based on fitness difference
                    edge_weight = exp(-(AllFitness(i) - AllFitness(j))^2);
                    connections(i,j) = edge_weight;
                    connections(j,i) = edge_weight;
                end
            end
        end
        last_update_iter = it;
    end
    
    %% 2. Enhanced Message Passing and Position Update
    % Calculate adaptive learning rate
    K = 1 - ((FEs)^(1/6)/(MaxFEs)^(1/6)); % From AO
    E = exp(-4*(FEs/MaxFEs));             % Exploitation weight
    
    for i = 1:N
        neighbors = find(connections(i,:));
        if ~isempty(neighbors)
            % Aggregate neighbor information (message passing)
            neighbor_influence = zeros(1, dim);
            total_weight = sum(connections(i,neighbors));
            
            for j = neighbors
                weight = connections(i,j) / total_weight;
                neighbor_influence = neighbor_influence + weight * (P(j,:) - P(i,:));
            end
            
            % Calculate exploration component (inspired by MGO)
            D_wind = w * (rand(1,dim)-0.5) * (1-FEs/MaxFEs);
            exploration = 0.1 * w * (rand(1,dim)-0.5) * (1-FEs/MaxFEs) * ...
                         (1 + 0.5*(1+tanh(beta/gamma))*(1-FEs/MaxFEs));
            
            % Combined update rule
            NewP = P(i,:) + ...
                  message_passing_strength * neighbor_influence + ... % GNN component
                  K * (Best_P - P(i,:)) + ...                       % AO component
                  E * (D_wind + exploration);                       % MGO component
            
            % Boundary handling
            NewP = min(max(NewP, lb), ub);
            
            % Evaluate and update
            new_fitness = fobj(NewP);
            FEs = FEs + 1;
            
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
    
    %% 3. Adaptive Diversity Management
    if mod(it, 10) == 0 && FEs < MaxFEs * 0.9
        [~, worst_idx] = max(AllFitness);
        % Hybrid reset strategy
        if rand() < 0.5
            P(worst_idx,:) = lb + (ub - lb) .* rand(1, dim);
        else
            P(worst_idx,:) = Best_P + 0.1 * (ub - lb) .* (rand(1, dim) - 0.5);
        end
        AllFitness(worst_idx) = fobj(P(worst_idx,:));
        FEs = FEs + 1;
    end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end