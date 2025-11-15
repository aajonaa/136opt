function [Best_P, Convergence_curve] = GNN_Fast2(N, MaxFEs, lb, ub, dim, fobj)
%% Initialize Parameters
Best_P = zeros(1, dim);
Best_fit = inf;
AllFitness = inf * ones(N, 1);
Convergence_curve = [];
FEs = 0;
it = 1;
% GNN parameters
update_prob = 1; % Increased from 0.3 for more active updates
k_neighbors = 5;   % Increased from 3 for better information flow
update_interval = 5;
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
    progressRatio = FEs / MaxFEs;
    
    %% 1. Graph Structure Update (maintaining GNN concept)
    if it - last_update_iter >= update_interval
        % Sort population by fitness for neighborhood structure
        [~, sorted_idx] = sort(AllFitness);
        connections = sparse(N, N);
        
        % Update graph structure using k-nearest neighbors
        for i = 1:N
            % Dynamic neighborhood size based on progress
            local_k = ceil(k_neighbors * (1 - 0.5 * progressRatio));
            start_idx = max(1, i - local_k);
            end_idx = min(N, i + local_k);
            potential_neighbors = sorted_idx(start_idx:end_idx);
            
            % Add connections with weighted edges
            for j = potential_neighbors'
                if i ~= j
                    % Edge weight based on fitness difference
                    fit_diff = abs(AllFitness(i) - AllFitness(j));
                    weight = 1 / (1 + fit_diff);
                    connections(i,j) = weight;
                    connections(j,i) = weight;
                end
            end
        end
        last_update_iter = it;
    end
    
    %% 2. Enhanced Message Passing and Position Update
    % Adaptive update probability
    current_update_prob = update_prob * (1 + progressRatio);
    update_mask = rand(N, 1) < current_update_prob;
    update_indices = find(update_mask);
    
    % Update selected nodes
    for i = update_indices'
        neighbors = find(connections(i,:));
        if ~isempty(neighbors)
            % Weight-based neighbor aggregation (message passing)
            weights = full(connections(i, neighbors));
            weights = weights / sum(weights);
            
            % Aggregate neighborhood information
            neighbor_info = zeros(1, dim);
            for n = 1:length(neighbors)
                neighbor_info = neighbor_info + weights(n) * P(neighbors(n),:);
            end
            
            % Enhanced solution generation inspired by high-performance version
            w1 = 1 - progressRatio^2;
            w2 = progressRatio;
            R = rand(1, dim);
            
            % Three candidate solutions using GNN-based information
            NewP1 = P(i,:) + R .* (w1 * (Best_P - P(i,:)) + w2 * (neighbor_info - P(i,:)));
            NewP2 = neighbor_info + R .* (w1 * (Best_P - P(i,:)) + w2 * (P(i,:) - neighbor_info));
            NewP3 = P(i,:) + R .* (neighbor_info - P(i,:)) + (1-R) .* (Best_P - P(i,:));
            
            % Boundary handling
            NewP1 = min(max(NewP1, lb), ub);
            NewP2 = min(max(NewP2, lb), ub);
            NewP3 = min(max(NewP3, lb), ub);
            
            % Evaluate candidates
            F_New1 = fobj(NewP1);
            F_New2 = fobj(NewP2);
            F_New3 = fobj(NewP3);
            FEs = FEs + 3;
            
            % Select best candidate
            [MinFit, idx] = min([F_New1, F_New2, F_New3]);
            if MinFit < AllFitness(i)
                switch idx
                    case 1
                        P(i,:) = NewP1;
                    case 2
                        P(i,:) = NewP2;
                    case 3
                        P(i,:) = NewP3;
                end
                AllFitness(i) = MinFit;
                
                if MinFit < Best_fit
                    Best_P = P(i,:);
                    Best_fit = MinFit;
                end
            end
        end
    end
    
    %% 3. Enhanced Diversity Through Graph Restructuring
    if mod(it, 20) == 0 && FEs < MaxFEs * 0.8
        [~, worst_idx] = max(AllFitness);
        
        % Reset worst solution with neighborhood influence
        neighbors = find(connections(worst_idx,:));
        if ~isempty(neighbors)
            % Combine best solution, random exploration, and neighborhood info
            neighbor_mean = mean(P(neighbors,:), 1);
            P(worst_idx,:) = 0.4 * Best_P + 0.3 * neighbor_mean + ...
                            0.3 * (lb + (ub - lb) .* rand(1, dim));
        else
            P(worst_idx,:) = 0.5 * Best_P + 0.5 * (lb + (ub - lb) .* rand(1, dim));
        end
        
        AllFitness(worst_idx) = fobj(P(worst_idx,:));
        FEs = FEs + 1;
    end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Positions = lb + (ub - lb) .* rand(SearchAgents_no, dim);
end