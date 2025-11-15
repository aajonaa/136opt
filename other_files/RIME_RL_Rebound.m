function [Best_rime, Convergence_curve] = RIME_RL_Rebound(N, MaxFEs, lb, ub, dim, fobj)
% 确保边界向量维度正确（核心修复）
if isscalar(lb)
    lb = lb * ones(1, dim);
else
    lb = lb(:)'; % 强制转换为行向量
    if length(lb) ~= dim
        error('Lower bound dimensions must match problem dimension');
    end
end
if isscalar(ub)
    ub = ub * ones(1, dim);
else
    ub = ub(:)'; % 强制转换为行向量
    if length(ub) ~= dim
        error('Upper bound dimensions must match problem dimension');
    end
end

% 初始化
Best_rime = zeros(1, dim);
Best_rime_rate = inf;
Rimepop = initialization(N, dim, ub, lb);
FEs = 0;
Convergence_curve = [];
Rime_rates = zeros(1, N);

% RL参数
learning_rate = 0.1;
discount_factor = 0.9;
epsilon = 1.0;
epsilon_decay = 0.999;
min_epsilon = 0.01;
num_states = 27;
num_actions = 3;
Q_table = zeros(num_states, num_actions);

% RIME参数
W = 5; a = 10; b = 0.5;
alpha_rebound = 0.5;
phi = 1.6;
history_size = 20;
best_rate_history = ones(1, history_size) * inf;

% 初始评估
for i = 1:N
    Rime_rates(i) = fobj(Rimepop(i, :));
    FEs = FEs + 1;
    if Rime_rates(i) < Best_rime_rate
        Best_rime_rate = Rime_rates(i);
        Best_rime = Rimepop(i, :);
    end
end
best_rate_history(1) = Best_rime_rate;

% 主循环
Time = 1;
while FEs < MaxFEs
    previous_best_rate = Best_rime_rate;
    
    % 状态获取（简化版）
    current_state = mod(floor(FEs/(MaxFEs/10)), num_states) + 1;
    
    % ε-贪婪动作选择
    if rand() < epsilon
        action = randi(num_actions);
    else
        [~, action] = max(Q_table(current_state, :));
    end
    
    % 计算RimeFactor
    progress = FEs/MaxFEs;
    RimeFactor = (rand - 0.5)*2*cos(pi*progress*10)*(1-round(progress*W)/W)/(1+exp(a*(progress-b)));
    E = sqrt(progress);
    
    % 安全归一化
    min_rate = min(Rime_rates);
    max_rate = max(Rime_rates);
    if abs(max_rate - min_rate) < eps
        normalized_rime_rates = zeros(size(Rime_rates));
    else
        normalized_rime_rates = (Rime_rates - min_rate)/(max_rate - min_rate);
    end
    
    newRimepop = Rimepop;
    
    % 种群更新（安全实现）
    for i = 1:N
        temp_pos = newRimepop(i, :);
        for j = 1:dim
            switch action
                case 1 % soft-rime
                    if rand() < E
                        temp_pos(j) = Best_rime(j) + RimeFactor*((ub(j)-lb(j))*rand + lb(j));
                    end
                case 2 % hard-rime
                    if rand() < normalized_rime_rates(i)
                        temp_pos(j) = Best_rime(j);
                    end
                case 3 % mix
                    if rand() < E
                        temp_pos(j) = Best_rime(j) + RimeFactor*((ub(j)-lb(j))*rand + lb(j));
                    end
                    if rand() < normalized_rime_rates(i)
                        temp_pos(j) = Best_rime(j);
                    end
            end
        end
        
        % 边界回弹（安全实现）
        for j = 1:dim
            if temp_pos(j) > ub(j)
                excess = temp_pos(j) - ub(j);
                temp_pos(j) = ub(j) - alpha_rebound * excess;
            elseif temp_pos(j) < lb(j)
                excess = lb(j) - temp_pos(j);
                temp_pos(j) = lb(j) + alpha_rebound * excess;
            end
        end
        
        % 评估新位置
        new_rate = fobj(temp_pos);
        FEs = FEs + 1;
        
        if new_rate < Rime_rates(i)
            Rime_rates(i) = new_rate;
            Rimepop(i, :) = temp_pos;
            if new_rate < Best_rime_rate
                Best_rime_rate = new_rate;
                Best_rime = temp_pos;
            end
        end
        
        if FEs >= MaxFEs
            break;
        end
    end
    
    % Q-learning更新
    if Best_rime_rate < previous_best_rate && previous_best_rate ~= inf
        reward = 10*(previous_best_rate - Best_rime_rate)/abs(previous_best_rate);
    else
        reward = -0.1;
    end
    
    next_state = mod(current_state, num_states) + 1;
    Q_table(current_state, action) = Q_table(current_state, action) + ...
        learning_rate*(reward + discount_factor*max(Q_table(next_state, :)) - Q_table(current_state, action));
    
    % ε衰减
    epsilon = max(min_epsilon, epsilon*epsilon_decay);
    
    % 记录历史
    best_rate_history = [best_rate_history(2:end), Best_rime_rate];
    Convergence_curve(Time) = Best_rime_rate;
    Time = Time + 1;
    
    if FEs >= MaxFEs
        break;
    end
end
end