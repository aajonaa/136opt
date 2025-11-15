%%      Optimization based on the smart behavior of plants with its engineering applications: %%
%%                               Ivy Algorithm (IVYA)                                         %%
%%                    Mojtaba Ghasemi, Mohsen Zare, Pavel Trojovský,                          %%
%%              Ravipudi Venkata Rao, Eva Trojovsk´,Venkatachalam Kandasamy                   %%
%%                           Knowledge-Based Systems (2024)                                   %%
%%                   DOI:https://doi.org/10.1016/j.knosys.2024.111850                         %%
%%                        https://www.optim-app.com/projects/ivya                             %%
%%             visit our website for more algorithms and their source code:                   %%
%%                             https://www.optim-app.com                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_pos,Convergence_curve]=IVYMod(N,MaxFEs,lb,ub,dim,fobj)

%% Initialization

X = zeros(N, dim);
AllFitness = inf * ones(N, 1);
newFitness = inf * ones(N, 1);
GV= zeros(N, dim);
newGV = zeros(N, dim);

FEs = 0;

for i = 1:N

    % Initialize Position
    %% Eq.(1)
    X(i, :) = unifrnd(lb, ub, 1, dim);
    %% Eq.(6)-upper condition
    GV(i, :)=(X(i, :) /(ub-lb));  

    % Evaluation
    AllFitness(i, 1) = fobj(X(i, :));
    FEs = FEs + 1;

end

% Initialize Best Cost History
BestFitness = zeros(MaxFEs, 1);

%% Ivy Main Loop
it = 1;
while FEs <= MaxFEs

    % Initialize new Population
    newX = [];

    for i = 1:N

        ii=i+1;
        if i==N
            ii=1;
        end

        beta_1=1+(rand/2); % beta value in line 8 of "Algorithm 1 (Fig.2 in paper)"

        if  AllFitness(i, 1)<beta_1*AllFitness(1, 1)
            %% Eq.(5)-(6)
            newSol =X(i, :) +abs(randn(1, dim)) .*(X(ii, :) - X(i, :))+ randn(1, dim) .*(GV(i, :));
        else
            %% Eq.(7)
            newSol =X(1, :) .*(rand+randn(1, dim).*(GV(i, :)));
        end
        
        %% Eq.(3)
        GV(i, :)=(GV(i, :)).*((rand^2)*(randn(1,dim)));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        newSol = max(newSol, lb);
        newSol = min(newSol, ub);
        %% Eq.(8)
        newGV(i, :)=(newSol /(ub-lb));

        %% Evaluate new population
        newFitness(i, 1) = fobj(newSol);
        FEs = FEs + 1;

        newX = [newX;
            newSol];
    end

    % Merge Populations
    X = [X;
        newX];
    Fitness = [AllFitness
        newFitness];

    % Sort Population
    [SortFitness, SortOrder]=sort(Fitness);
    X = X(SortOrder, :);
    AllFitness = SortFitness(1:N, 1);

    % Competitive Exclusion (Delete Extra Members)
    if numel(X)>N
        X = X(1:N, :);
    end

    % Store Best Solution Ever Found
    best_pos = X(1);
    BestFitness(FEs) = AllFitness(1);
    Convergence_curve(it) =AllFitness(1);
    it = it + 1;
end


%% Results
bestFitness=Fitness(1);
best_pos=X(1, :);
end