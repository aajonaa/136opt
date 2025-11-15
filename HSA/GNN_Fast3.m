function [Best_P, Convergence_curve] = GNN_Fast3(N, MaxFEs, lb, ub, dim, fobj)
%% Initialize Parameters
Best_P = zeros(1, dim);
Best_fit = inf;
AllFitness = inf * ones(N, 1);
Convergence_curve = [];
FEs = 0;
it = 1;

% GNN parameters
message_passing_strength = 0.3;
attention_weight = 0.5;
graph_sparsity = 4;
update_interval = 3;
theta = 0.5; % Connection threshold from high-performance version

% Initialize graph structure
connections = sparse(N, N);
similarity = zeros(N, N);
NoUpdateCount = zeros(1, N);
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
    
    %% 1. Graph Structure Update with Attention
    if it - last_update_iter >= update_interval
        connections = sparse(N, N);
        [SortFit, FitRank] = sort(AllFitness);
        
        for i = 1:N
            for j = i+1:N
                % Calculate similarity (from high-performance version)
                NormDistIJ = sqrt(sum(((P(i,:) - P(j,:)) ./ (ub - lb)) .^ 2)) / sqrt(dim);
                DiffFitIJ = abs(AllFitness(i) - AllFitness(j)) / (max(AllFitness) - min(AllFitness));
                
                % Calculate similarity with attention weight
                similarity(i,j) = 1 - log(1 + ((1-progressRatio) * NormDistIJ + progressRatio * DiffFitIJ) * (exp(1) - 1));
                similarity(j,i) = similarity(i,j);
                
                % Establish connections based on similarity threshold
                if similarity(i,j) >= theta
                    connections(i,j) = 1;
                    connections(j,i) = 1;
                end
            end
        end
        last_update_iter = it;
    end
    
    %% 2. Message Passing and Position Update
    % Adaptive adjustment factor from high-performance version
    AAF = 0.5 * 2^(1 - progressRatio);
    
    for i = 1:N
        neighbors = find(connections(i,:));
        if ~isempty(neighbors)
            % Calculate local mean (from high-performance version)
            LocalMean = mean(P(neighbors,:), 1);
            
            % Generate candidate solutions using Levy flight
            R = Levy(1,dim) .* (rand(1,dim) > 0.5);
            w1 = 1 - progressRatio^2;
            
            % Solution generation strategies from high-performance version
            Node1 = P(i,:) + AAF * R .* (w1 * (Best_P - LocalMean) + (1-w1) * (LocalMean - P(i,:)));
            Node2 = Best_P + AAF * R .* (w1 * (LocalMean - P(i,:)) + (1-w1) * (P(i,:) - Best_P));
            Node3 = LocalMean + AAF * R .* (w1 * (P(i,:) - Best_P) + (1-w1) * (Best_P - LocalMean));
            
            % Boundary handling
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
    
    %% 3. Population Diversity Management (from high-performance version)
    % if mod(it, 10) == 0
        [~, pruneIdx] = maxk(NoUpdateCount, 1 + floor(rand * N * 0.1));
        for IDX = pruneIdx
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
    % end
    
    Convergence_curve(it) = Best_fit;
    it = it + 1;
end
end

% Levy flight (from high-performance version)
function o = Levy(n,d)
    beta = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(n,d) * sigma;
    v = randn(n,d);
    step = u./abs(v).^(1/beta);
    o = step;
end

function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Positions = lb + (ub - lb) .* rand(SearchAgents_no, dim);
end