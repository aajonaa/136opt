
% =============================================================================================================================================================================
% Differentiated Creative Search (DCS)
% CITATION:
% Duankhan P., Sunat K., Chiewchanwattana S., and Nasa-ngium P. "The Differentiated Creative Search (DCS): Leveraging Differentiated Knowledge-acquisition and Creative Realism
% to Address Complex Optimization Problems". (Accepted for publication in Expert Systems with Applications)
% =============================================================================================================================================================================
function [best_pos,Convergence_curve] = DCS(N,MaxFEs,lb,ub,dim,fobj)
    lb = lb * ones(1, dim);
    ub = ub * ones(1, dim);
    % Parameters
    newX = zeros(N,dim);
    new_pos = zeros(N,dim);
    Convergence_curve = [];
    eta_qKR = zeros(1,N); %%%%%
    newAllFitness = zeros(N,1);
    % Golden ratio
    golden_ratio = 2/(1 + sqrt(5)); %%%%%
    % High-performing individuals
    ngS = max(6,round(N * (golden_ratio/3))); %%%%%

    % Initialize the population
    X = zeros(N,dim);
    for i = 1:N
        X(i,:) = lb + rand(1,dim) .* (ub - lb);
    end

    % Initialize fitness values
    AllFitness = zeros(N,1);
    for i = 1:N
        AllFitness(i,1) = fobj(X(i, :));
    end
    % Generation
    FEs = 0;
    it = 1;
    pc = 0.5; %%%%%
    % Best solution
    bestFitness = min(AllFitness);
    % Ranking-based self-improvement
    phi_qKR = 0.25 + 0.55 * ((0 + ((1:N)/N)) .^ 0.5); %%%%%
    while FEs < MaxFEs
        % Sort population by fitness values
        [X, AllFitness, ~] = PopSort(X,AllFitness);
        % Reset
        bestInd = 1;
        % Compute social impact factor
        lamda_t = 0.1 + (0.518 * ((1-(FEs/MaxFEs)^0.5))); %%%%%
        for i = 1:N
            % Compute differentiated knowledge-acquisition rate
            eta_qKR(i) = (round(rand * phi_qKR(i)) + (rand <= phi_qKR(i)))/2; %%%%%
            jrand = floor(dim * rand + 1); %%%%%
            newX(i,:) = X(i,:);
            if i == N && rand < pc
                % Low-performing
                newX(i,:) = lb + rand * (ub - lb);
            elseif i <= ngS %%%%%
                % High-performing
                while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= bestInd, break, end, end
                for d = 1:dim
                    if rand <= eta_qKR(i) || d == jrand
                        newX(i,d) = X(r1,d) + LnF3(golden_ratio,0.05,1,1);
                    end
                end
            else
                % Average-performing
                while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= bestInd, break, end, end
                while true, r2 = ngS + round((N - ngS) * rand + 0.5); if r2 ~= i && r2 ~= bestInd && r2 ~= r1, break, end, end
                % Compute learning ability
                omega_it = rand; %%%%%
                for d = 1:dim
                    if rand <= eta_qKR(i) || d == jrand
                        newX(i,d) = X(bestInd,d) + ((X(r2,d) - X(i,d)) * lamda_t) + ((X(r1,d) - X(i,d)) * omega_it);
                    end
                end
            end
            % Boundary
            newX(i,:) = boundConstraint(newX(i,:),X(i,:),[lb; ub]);
            new_pos(i, :) = newX(i, :);
            newAllFitness(i,1) = fobj(newX(i, :));
            FEs = FEs + 1;
            if newAllFitness(i,1) <= AllFitness(i,1)
                X(i,:) = new_pos(i,:);
                AllFitness(i,1) = newAllFitness(i,1);
                if newAllFitness(i,1) < bestFitness
                    bestFitness = newAllFitness(i,1);
                    bestInd = i;
                end
            end
        end
        best_pos = X(bestInd,:);
        best_cost = bestFitness;
        Convergence_curve(it) = bestFitness;
        it = it + 1;
    end
end


function [sorted_population, sorted_fitness, sorted_index] = PopSort(input_pop,input_fitness)
    [sorted_fitness, sorted_index] = sort(input_fitness,1,'ascend');
    sorted_population = input_pop(sorted_index,:);
end


function Y = LnF3(alpha, sigma, m, n)
    Z = laplacernd(m, n);
    Z = sign(rand(m,n)-0.5) .* Z;
    U = rand(m, n);
    R = sin(0.5*pi*alpha) .* tan(0.5*pi*(1-alpha*U)) - cos(0.5*pi*alpha);
    Y = sigma * Z .* (R) .^ (1/alpha);
end


function x = laplacernd(m, n)
    u1 = rand(m, n);
    u2 = rand(m, n);
    x = log(u1./u2);
end


function vi = boundConstraint(vi, pop, lu)
    % if the boundary constraint is violated, set the value to be the middle
    % of the previous value and the bound
    %
    % Version: 1.1   Date: 11/20/2007
    % Written by Jingqiao Zhang, jingqiao@gmail.com
    [NP, ~] = size(pop);  % the population size and the problem's dimension
    % check the lower bound
    xl = repmat(lu(1, :), NP, 1);
    pos = vi < xl;
    vi(pos) = (pop(pos) + xl(pos)) / 2;
    % check the upper bound
    xu = repmat(lu(2, :), NP, 1);
    pos = vi > xu;
    vi(pos) = (pop(pos) + xu(pos)) / 2;
end