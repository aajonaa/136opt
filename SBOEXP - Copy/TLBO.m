function [best_pos,Convergence_curve] = TLBO(N, Max_FEs, lb, ub, dim, fobj)

    %% Initialization 
    FEs = 0;
    
    % Initialize Population Array
    X = initialization(N, dim, ub, lb);

    % Initialize Best Solution
    bestFitness = inf;
    AllFitness = inf * ones(N, 1);
    X2_Fitness = inf * ones(N, 1);
    X3_Fitness = inf * ones(N, 1);

    % Initialize Population Members
    for i = 1:N
        AllFitness(i, 1) = fobj(X(i, :));
        FEs = FEs + 1;

        if AllFitness(i, 1) < bestFitness
            bestFitness = AllFitness(i, 1);
            best_pos = X(i, :);
        end
    end

    Convergence_curve = [];
    it = 1;

    %% TLBO Main Loop

    while FEs < Max_FEs

        % Calculate Population Mean
        Mean = 0;
        for i = 1:N
            Mean = Mean + X(i, :);
        end
        Mean = Mean/N;

        % Select Teacher
        Teacher = best_pos;

        % Teacher Phase

        % Create Empty Solution

        for i = 1:N

            % Teaching Factor
            TF = randi([1 2]);

            % Teaching (moving towards teacher)
            temp = X(i, :) + rand() .* (Teacher - TF*Mean);

            % Clipping
            Flag4ub=temp>ub;
            Flag4lb=temp<lb;
            temp=(temp.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

            % Evaluation
            tempFitness = fobj(temp);
            FEs = FEs + 1;

            % Comparision
            if tempFitness<AllFitness(i, 1)
                AllFitness(i ,1) = tempFitness;
                X(i, :) = temp;
                if AllFitness(i, 1) < bestFitness
                    bestFitness = AllFitness(i, 1);
                    best_pos = X(i, :);
                end
            end
        end

        % Learner Phase
        for i = 1:N

            A = 1:N;
            A(i) = [];
            j = A(randi(N-1));

            Step = X(i, :) - X(j, :);
            if AllFitness(j, 1) < AllFitness(i, 1)
                Step = -Step;
            end

            % Teaching (moving towards teacher)
            temp = X(i, :) + rand() .* Step;

            % Clipping
            Flag4ub=temp>ub;
            Flag4lb=temp<lb;
            temp=(temp.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

            % Evaluation
            tempFitness = fobj(temp);
            FEs = FEs + 1;

            % Comparision
            if tempFitness<AllFitness(i, 1)
                AllFitness(i, 1) = tempFitness;
                X(i, :) = temp;
                if AllFitness(i, 1) < bestFitness
                    bestFitness = AllFitness(i, 1);
                    best_pos = X(i, :);
                end
            end
        end

        % Store Record for Current Iteration
        Convergence_curve(it) = bestFitness;
        it = it + 1;

    end
    toc
end
