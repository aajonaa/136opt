%% Jona 2024-11-4 2024-12-2 2024-12-3 2024-12-12
% 2024-12-13: only if j == Best_id can not make all the nodes connect to
% the best nodes, so I changed it to if i == Best_id || j == Best_id
% 2024-12-14: Check the network connect status
% 2024-12-17: There have a error in the node cooperation search, the second
% formulation, I modify it
% 2024-12-17: Drop the connection when drop the nodes in the second phase
% 2024-12-18 Only take connected neighbor as consideration 
% 2024-12-18 No need to connected to the best node with each node
function [Best_P, Convergence_curve] = GNN(N,MaxFEs,lb,ub,dim,fobj)
%% 初始化参数
Best_P = zeros(1,dim);
Best_fit = inf;
AllFitness = inf*ones(N,1);
Convergence_curve=[];
theta = 0.5; % 连接阈值
maxNormDist = sqrt(dim); % 归一化距离的最大可能值
similarity = zeros(N, N); % 存储相似度
connections = zeros(N, N); % 保存连接的状态和相似度
NoUpdateCount = zeros(1,N); % 当前位置未更新的次数
it=1;
FEs=0;
%% 步骤1：节点（种群）初始化
P = initialization(N,dim,ub,lb);
for i=1:N
    FEs = FEs+1;
    AllFitness(i) = fobj(P(i,:));
    if AllFitness(i) < Best_fit
        Best_P = P(i,:);
        Best_fit = AllFitness(i);
    end
end
%%
while FEs <= MaxFEs
    % 计算当前评估的比例
    progressRatio = FEs / MaxFEs; 

    %% a. 节点连接构建  --->  其实是，通过位置差异或者适应度差异，分成动态子种群
    [SortFit,FitRank] = sort(AllFitness);
    Best_id = FitRank(1);
    Wrost_fit = SortFit(end);
    for i = 1:N
        for j = i+1:N
            %% Calculate similarity between node i and node j
            similarity = CalculateSimilarity(similarity, P, i, j, AllFitness, lb, ub, FEs, MaxFEs);

            %% Connect nodes between node i and node j
            connections = NodesConnect(connections, similarity, i, j, Best_id, theta);
        end
    end

    %% b. 节点重塑（Node Reshaping）--->增加多样性
    % Fitness-Weighted Centrality 适应度加权中心性 -- jona add
    FWC = calculateFWC(P, AllFitness, similarity, connections, lb, ub); % FWC的值越大，代表节点越核心
    [~, FWCrank] = sort(FWC, 'descend');
    % 定义剪枝指标Ci，Pi质量越好时，Ci越小
    Ci = zeros(N, 1);
    for i = 1:N
        Ci(i) = 0.6 * FitRank(i)/N + 0.3 * NoUpdateCount(i)./(max(NoUpdateCount) + eps) + 0.1 * (FWCrank(i)/ N);
    end
    % 选择剪枝的节点
    [~, pruneIndices] = maxk(Ci, 1 + floor(rand*N*0.2));  
    % 进行节点剪枝与替换
    for idx = pruneIndices'
        P(idx, :) = lb + (ub - lb) .* rand(1, dim); 

        % Drop the connection alos 12-17
        connections(idx, :) = 0;

        AllFitness(idx) = fobj(P(idx, :)); 
        FEs = FEs + 1;
        NoUpdateCount(idx) = 0; %
        if  AllFitness(idx) < Best_fit
            Best_P = P(idx, :);
            Best_fit = AllFitness(idx);
        end
        if AllFitness(idx) > Wrost_fit
            Wrost_fit =  AllFitness(idx);
        end
    end

    %% c.节点集探索（Node-set exploration）
    % Aaptive Adjustment Factor 自适应调节因子 jona-add
    AAF = 0.5 * 2^(1 - (FEs / MaxFEs));
    NewP = P;
    for i = 1:N
        if i < N
            %% c.1 节点邻居竞争 --->局部增强
            for j = i+1:N
                if connections(i, j) > 0 % 检查是否存在超边连接解i和解j
                    direction = sign(AllFitness(i)-AllFitness(j)).*ones(1,dim);
                    % 基于节点权重和随机扰动生成新解
                    NewP(i,:) = P(i,:) + (P(j,:) - P(i,:)) .* direction * AAF;
                    % 保证新解符合上下界限制
                    NewP(i,:) = min(max(NewP(i,:), lb), ub);
                    % 评估新解的适应度
                    f_pnew = fobj(NewP(i,:));
                    FEs = FEs + 1;
                    % 找出与当前超边相关的解中质量更差的解
                    index = [i,j];
                    [~, idx] = sort([AllFitness(i), AllFitness(j)]);
                    worst_idx = index(idx(2));
                    % 如果新点的适应度优于最差解的适应度，则进行替换，并更新超边权重
                    if f_pnew < AllFitness(worst_idx)
                        NoUpdateCount(worst_idx) = 0;
                        P(worst_idx, :) = NewP(i,:);
                        AllFitness(worst_idx) = f_pnew;
                        if f_pnew < Best_fit
                            Best_P = NewP(i,:);
                            Best_fit = f_pnew;
                        end
                    else
                        NoUpdateCount(worst_idx) = NoUpdateCount(worst_idx) + 1;
                    end
                end
            end
        end

        %% c.2 节点协同搜索--->全局探索
        LocalPosition = zeros(1, dim); 
        for j = 1:N
            if connections(i, j) > 0 % 考虑基于Pi超边上有连接的解
                LocalPosition = LocalPosition + P(j,:);
            end
        end
        LocalMean = LocalPosition./sum(connections(i, :));
        % 通过调节当前节点位置、局部平均位置与全局最优解之间的相互作用，实现解空间的有效探索与利用
        w1 = 1 - progressRatio^2;
        R = Levy(1,dim).*(rand(1,dim)>0.5);
        Node1 = P(i,:) + AAF * R .* (w1 * (Best_P - LocalMean) + (1-w1) * (LocalMean - P(i,:)));
        Node2 = Best_P + AAF * R .* (w1 * (LocalMean - P(i,:)) + (1-w1) * (P(i,:) - Best_P));
        Node3 = LocalMean +AAF * R .* (w1 * (P(i,:) - Best_P) + (1-w1) * (Best_P - LocalMean));
        % 保证新解符合上下界限制
        Node1 = min(max(Node1, lb), ub);
        Node2 = min(max(Node2, lb), ub);
        Node3 = min(max(Node3, lb), ub);
        % 计算适应度值
        F_Node1 = fobj(Node1);
        F_Node2 = fobj(Node2);
        F_Node3 = fobj(Node3);
        FEs=FEs+3;
        Fitness_comb=[F_Node1,F_Node2,F_Node3,AllFitness(i)];
        [Min_fitness_comb,m]=min(Fitness_comb);
        AllFitness(i) = Min_fitness_comb;
        if m==1   
            P(i,:)=Node1;
            NoUpdateCount(i) = 0;
        elseif m==2
            P(i,:)=Node2; 
            NoUpdateCount(i) = 0;
        elseif m==3
            P(i,:)=Node3; 
            NoUpdateCount(i) = 0;
        else
            P(i,:)=P(i,:); 
            NoUpdateCount(i) = NoUpdateCount(i) + 1;
        end
        if AllFitness(i) < Best_fit
           Best_P = P(i,:);
           Best_fit = AllFitness(i);
        end
    end
    %%
    Convergence_curve(it)=Best_fit;
    it=it+1;
end
end

% Similaritry calculation
function similarity = CalculateSimilarity(similarity, P, i, j, AllFitness, lb, ub, FEs, MaxFEs)
    % Normalized distance of node i and node j
    NormDistIJ = sqrt(sum(((P(i, :) - P(j, :)) ./ (ub - lb)) .^ 2)) / sqrt(size(P, 2));
    % NormDistIJ = norm((P(i, :) - P(j, :)) ./ (ub - lb)) / sqrt(size(P, 2));
    % Fitness difference between node i and node j
    DiffFitIJ = abs(AllFitness(i) - AllFitness(j)) / (max(AllFitness) - min(AllFitness));
    % Similarity of node i and node j (Distribute different
    % weight to the NormDistIJ and DiffFitIJ: self adaptive
    % weight)
    similarity(i, j) = 1 - log(1 + ((1-FEs/MaxFEs) * NormDistIJ + FEs/MaxFEs * DiffFitIJ) * (exp(1) - 1));
    % similarity(i, j) = 1 - (1-FEs/MaxFEs) * NormDistIJ - FEs/MaxFEs * DiffFitIJ;
    similarity(j ,i) = similarity(i, j);
    similarity(i, i) = 0.1;
end

%% Nodes connection
function connections = NodesConnect(connections, similarity, i, j, Best_id, theta)
    % if i == Best_id || j == Best_id
    %     connections(i, j) = 1;
    %     connections(j, i) = 1;
    % else
        if similarity(i, j) >= theta
            connections(i, j) = 1;
            connections(j, i) = 1;
        end
    % end
end

function FWC = calculateFWC(P, AllFitness, similarity, connections, lb, ub)
    N = size(P, 1); % 种群大小
    FWC = zeros(N, 1); % 初始化每个节点的FWC值
    epsilon = 1e-4; % 避免除以零
    for i = 1:N
        % 计算节点i与其所有邻居的加权距离之和
        weightedDistSum = 0;
        for j = 1:N
            if connections(i, j) > 0 % 仅考虑存在相似度的邻居 % Only take connected neighbor as consideration 2024-12-18
                d_ij = norm(P(i,:) - P(j,:)) / norm(ub - lb); % 归一化欧式距离
                weightedDistSum = weightedDistSum + similarity(i, j) / (d_ij + epsilon);
            end
        end
        FWC(i) = 1 / (AllFitness(i) + epsilon) + weightedDistSum;
    end
end

function o = Levy(n,d)
beta=1.5;
sigma=(gamma(1+beta).*sin(pi*beta/2)./(gamma((1+beta)/2).*beta.*2.^((beta-1)/2))).^(1/beta);
u=randn(n,d)*sigma;
v=randn(n,d);
step=u./abs(v).^(1/beta);
o=step;
end