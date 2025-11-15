% improved Dual Dispersal mechanisms for Status-based Optimization (DDSBO)
% To prevent SBO from premature convergence in the early stage and improve its ability to escape from local search, 
% this study combines the superior dispersal mechanisms of the animated oat optimization algorithm (AAO) and moss 
% growth optimization (MGO), and proposes the improved Dual Dispersal mechanisms for Status-based Optimization (DDSBO).
function [best_pos, Convergence_curve] = DDSBO(N, MaxFEs, lb, ub, dim, fobj)

%% INITIALIZATION
FEs = 0;
bestFitness = inf; % Change to -inf for maximization problems
best_pos = zeros(1, dim);

Convergence_curve = [];
iter = 1;

if length(lb) == 1
    lb = lb * ones(1, dim);
end
if length(ub) == 1
    ub = ub * ones(1, dim);
end


% Initialize populations
current_X = initialization(N, dim, ub, lb);
localElite_X = initialization(N, dim, ub, lb);

% Initialize fitness values
current_Fitness = inf * ones(N, 1);
localElite_Fitness = inf * ones(N, 1);

% Calculate initial fitness and determine initial best solutions (local and global)
for i = 1:N
    current_Fitness(i, 1) = fobj(current_X(i, :));
    FEs = FEs + 1;

    temp_localElite_Fitness = fobj(localElite_X(i ,:));
    FEs = FEs + 1;

    % Determine the initial local elite
    if current_Fitness(i, 1) < temp_localElite_Fitness
        localElite_X(i, :) = current_X(i, :);
        localElite_Fitness(i, 1) = current_Fitness(i, 1);
    else
        localElite_Fitness(i, 1) = temp_localElite_Fitness;
    end

    % Update the initial global best
    if localElite_Fitness(i, 1) < bestFitness
        bestFitness = localElite_Fitness(i, 1);
        best_pos = localElite_X(i, :);
    end
end

%% Sort the local elite fitness for the first Roulette Wheel Selection
[sorted_localElite_Fitness, ~] = sort(localElite_Fitness);

%% Social success flag and social fitness initialization
flag = ones(N, 1);
social_Fitness = inf * ones(N, 1);

m = 0.5 * rand(1, N) / dim;
L = N / dim * rand(1, N);
e = m;
g = 9.8/dim;
x = 3 * rand(1, N) / dim;

%% MAIN LOOP
while FEs < MaxFEs


    %% Select an individual from the localElite population based on Roulette selection
    Roulette_index = RouletteWheelSelection(1./(sorted_localElite_Fitness + eps));
    if Roulette_index == -1
        Roulette_index = 1; % Default to the best if selection fails
    end

    

    %% Update the current population
    for i = 1:N
        w1 = randn;
        w2 = randn;
        w3 = tanh((sqrt(abs(MaxFEs - randn * FEs))/i)^(FEs/MaxFEs));
        w4 = unifrnd(-w3, w3);
        if rand < w3
            for j = 1:dim
                current_X(i, j) = (1 - w1 - w2) * current_X(i, j) + w1 * localElite_X(Roulette_index, j) + w2 * best_pos(j);
            end
        else
            for j = 1:dim
                current_X(i, j) = w4 * ((1 - w1 - w2) * current_X(i, j) + w1 * localElite_X(Roulette_index, j) + w2 * best_pos(j));
            end
        end
        % if FEs/MaxFEs * rand > 0.1
            % current_X(i,:) = current_X(i,:) + step .* D_wind;
        % end
        m = zeros(1, dim);
        u = randperm(dim);
        m(u(1:ceil(rand * dim))) = 1;
        c = (1 -  FEs / MaxFEs)^3;
        P = levy(N,dim,1.5);
        theta = pi * rand(1,N);

        if rand * FEs/MaxFEs > 0.4
            W = c/pi *(2 * rand(1,dim) - 1) .* ub  ;
            if mod(i,N/10) == 0
                current_X(i,:) = mean(current_X) + W;
            elseif mod(i,N/10) == 1
                current_X(i,:) = best_pos + W ;
            else
                current_X(i,:) = current_X(i,:) + W ;
            end

            % if rand > 0.5
                % A = ub - abs(ub * FEs * sin(2 * pi * rand) / MaxFEs);
                % R = (m(i) * e(i) + L(i) ^2) /dim * unifrnd( -A , A, 1, dim);
                % current_X(i,:) = best_pos + R  + c * P(i,:) .* best_pos;
            % else
            %     k = 0.5 + 0.5 * rand;
            %     B = ub - abs(ub * FEs * cos(2 * pi * rand) / MaxFEs);
            %     alpha = 1 / pi * exp((randi([0,FEs]) / MaxFEs));
            %     J = 2 * k * x(i)^2 * sin (2 * theta(i)) / m(i) / g * (1 - alpha) /dim  * unifrnd( -B , B, 1, dim);
            %     current_X(i,:) = best_pos + J + c * P(i,:) .* best_pos;
            % end
        end
    end

    calPositions = current_X;
    div_num = randperm(dim);
    divide_num = dim/4;
    m_w = 2;
    %Divide the population and select the regions with more individuals based on the best
    for j=1:max(divide_num,1)
        th = best_pos(div_num(j));
        index = calPositions(:,div_num(j)) > th;
        if sum(index) < size(calPositions, 1)/2 %choose the side of the majority
            index = ~index;
        end
        calPositions = calPositions(index,:);
    end
    D = best_pos - calPositions; %Compute the distance from individuals to the best
    D_wind = sum(D, 1)/size(calPositions, 1); %Calculate the mean of all distances
    step = m_w * (rand(size(D_wind))-0.5) * (1-FEs/MaxFEs);

    for i = 1:N
        current_X(i,:) = current_X(i,:) + step .* D_wind;
    end





    %% Boundary control
    current_X = BoundaryControl(current_X, lb, ub);

    %% Upward social strategy
    social_X = current_X;

    % One-dimension source exchange for socially successful individuals
    for i = 1:N
        if flag(i) == 1
            social_X1_val = localElite_X(i, randi(dim));
            social_X2_val = best_pos(randi(dim));
            social_X(i, randi(dim)) = (social_X1_val + social_X2_val) / 2;
        end
    end

    % Multi-dimension source exchange for socially failed individuals

    for i = 1:N
        if flag(i) == 0
            for j = 1:dim
                if m(j)
                    social_X(i, j) = localElite_X(i, j);
                end
            end
        end
    end


    %% Greedy selection and current population update
    for i = 1:N
        % Evaluate new and social positions
        if FEs < MaxFEs
            current_Fitness(i, 1) = fobj(current_X(i, :));
            FEs = FEs + 1;
        end
        if FEs < MaxFEs
            social_Fitness(i, 1) = fobj(social_X(i, :));
            FEs = FEs + 1;
        end

        % Greedy selection
        if social_Fitness(i, 1) < current_Fitness(i, 1)
            % Social success: update position and set flag
            flag(i, 1) = 1;
            current_X(i, :) = social_X(i, :);
            current_Fitness(i, 1) = social_Fitness(i, 1);
        else
            % Social fail: keep current position and set flag
            flag(i, 1) = 0;
        end
    end

    %% Update local elite population
    for i = 1:N
        if current_Fitness(i, 1) < localElite_Fitness(i, 1)
            localElite_Fitness(i, 1) = current_Fitness(i, 1);
            localElite_X(i, :) = current_X(i, :);
        end
    end

    %% Sort local elite fitness and update the global best position
    [sorted_localElite_Fitness, indBest] = sort(localElite_Fitness);
    if sorted_localElite_Fitness(1) < bestFitness
        bestFitness = sorted_localElite_Fitness(1);
        best_pos = localElite_X(indBest(1), :);
    end

    %% Record the best fitness for the convergence curve
    Convergence_curve(iter) = bestFitness;
    iter = iter + 1;
end
end

%% Helper Functions

% Enforce boundary constraints on agent positions
function X = BoundaryControl(X, low, up)
[N, dim] = size(X);

if isscalar(low)
    low = repmat(low, 1, dim);
end
if isscalar(up)
    up = repmat(up, 1, dim);
end

for i = 1:N
    for j = 1:dim
        k = rand < rand; % 50% chance for either clipping or re-initializing

        if X(i,j) < low(j)
            if k
                X(i,j) = low(j); % Clipping to the boundary
            else
                X(i,j) = rand * (up(j) - low(j)) + low(j); % Re-initializing
            end
        end

        if X(i,j) > up(j)
            if k
                X(i,j) = up(j); % Clipping to the boundary
            else
                X(i,j) = rand * (up(j) - low(j)) + low(j); % Re-initializing
            end
        end
    end
end
end

% Roulette wheel selection mechanism
function choice = RouletteWheelSelection(weights)
% Normalize weights to handle negative values if any, although 1/fitness should be positive
if any(weights<0)
    weights = weights + min(weights);
end
accumulation = cumsum(weights);
if accumulation(end) == 0
    choice = -1;
    return;
end
p = rand() * accumulation(end);
chosen_index = -1;
for index = 1:length(accumulation)
    if (accumulation(index) > p)
        chosen_index = index;
        break;
    end
end
choice = chosen_index;
end

function [z] = levy(n,m,beta)

num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator

den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator

sigma_u = (num/den)^(1/beta);% Standard deviation

u = random('Normal',0,sigma_u,n,m);

v = random('Normal',0,1,n,m);

z =u./(abs(v).^(1/beta));


end