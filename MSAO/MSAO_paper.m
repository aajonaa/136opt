%% 2024-7-1
function [best_pos, Convergence_curve] = MSAO(N, Max_FEs, lb, ub, dim, fobj)
    % function [bestFitness, Convergence_curve]=AO(N,Max_FEs,lb,ub,dim,fobj) % For myself

    FEs = 0;
    it = 1;
    Fitn = inf * ones(N, 1);
    newFitAll = inf * ones(N, 1);
    FitTrack = inf * ones(N, 3);
    FitAll = inf * ones(N, 1);
    Convergence_curve = [];

    X = initialization(N, dim, ub, lb);

    % Initial fitness calculation and fitTrack initialization
    for i = 1:N
        FitAll(i) = fobj(X(i, :));
    end
    FEs = FEs + N;
    FitTrack = repmat(FitAll, 1, 3);

    [bestFitness, minfit_idx] = min(FitAll);
    best_position = X(minfit_idx, :);
    newX = zeros(N, dim);

    while FEs <= Max_FEs
        P1 = 1 - ((FEs)^(1/3) / (Max_FEs)^(1/3));
        d = 1 * exp(-4 * (FEs / Max_FEs));
        TDR = P1;

        FitTrackNorm = (FitTrack - min(FitTrack)) ./ (max(FitTrack) - min(FitTrack));
        alpha = (FitTrackNorm(:, 2) - FitTrackNorm(:, 1)) ./ (sum(FitTrackNorm(:, 2) - FitTrackNorm(:, 1)) + eps) + 0.2;
        factor = 1/2;

        for i = 1:N
            Fitn(i) = (FitAll(i) - min(FitAll)) / (max(FitAll) - min(FitAll));
            for j = 1:dim
                if rand < P1
                    [~, index] = sort(FitAll);
                    NumOfAlpha = floor(N / 3);
                    tempIdx = index(1:NumOfAlpha);
                    choose = randperm(NumOfAlpha);
                    newX(i, j) = X(i, j) - d * (X(randi(N), j) - X(i, j)) + d * (X(tempIdx(choose(1)), j) - X(i, j));
                else
                    theta = pi/2 * min(Fitn) / Fitn(i);
                    newX(i, j) = X(i, j) + sin(theta) * best_position(j);
                end
            end

            % Boundary control
            newX(i, :) = Transborder_reset(newX(i, :), lb, ub, dim, best_position);

            % Update fitness and apply greedy selection
            tempFit = fobj(newX(i, :));
            FEs = FEs + 1;
            if tempFit < FitAll(i)
                X(i, :) = newX(i, :);
                FitAll(i) = tempFit;
            end

            % Two-stage mutation and boundary control
            if FEs < factor * Max_FEs
                newX(i, :) = Mutation(newX(i, :), X(i, :), best_position, dim);
            else
                newX(i, :) = best_position + alpha(i) * TDR * randsrc(1, dim) .* ((ub - lb) .* rand(1, dim) + lb);
            end
            newX(i, :) = Transborder_reset(newX(i, :), lb, ub, dim, best_position);
            newFitAll(i) = fobj(newX(i, :));
            FEs = FEs + 1;
        end

        % Update fitness track and apply greedy selection
        FitTrack = [newFitAll, FitTrack(:, 1:2)];
        for i = 1:N
            if newFitAll(i) < FitAll(i)
                X(i, :) = newX(i, :);
                FitAll(i) = newFitAll(i);
            end
        end

        % Update best position
        [fit_min, minfit_idx] = min(FitAll);
        if fit_min < bestFitness
            best_position = X(minfit_idx, :);
            bestFitness = fit_min;
        end

        Convergence_curve(it) = bestFitness;
        it = it + 1;
    end

    best_pos = best_position;
end

function z = Mutation(z, x, b, dim)
    for j = 1:dim
        if rand < 0.05
            z(j) = x(j);
        end
        if rand < 0.2
            z(j) = b(j);
        end
    end
end

function z = Transborder_reset(z, lb, ub, dim, best)
    for j = 1:dim
        if z(j) > ub || z(j) < lb || isnan(z(j))
            z(j) = best(j);
        end
    end
end
