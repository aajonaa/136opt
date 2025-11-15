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
function [best_pos,Convergence_curve]=IVYMod2(N,Max_FEs,lb,ub,dim,fobj)

%% Initialization

% Empty Plant Structure
% empty_plant.Position = [];
% empty_plant.Cost = [];
% empty_plant.GV= [];
Cost = inf * ones(N, 1);
newCost = inf * ones(N, 1);
GV= zeros(N, dim);
newGV= zeros(N, dim);
pop = zeros(N, dim);    % Initial Population Array
newpop = zeros(N, dim);    % Initial Population Array

FEs = 0;

for i = 1:N

    % Initialize Position
    %% Eq.(1)
    pop(i, :) = unifrnd(lb, ub, 1, dim);
    %% Eq.(6)-upper condition
    GV(i, :)=(pop(i, :) /(ub-lb));  

    % Evaluation
    Cost(i) = fobj(pop(i, :));
    FEs = FEs + 1;

end

[~, index] = sort(Cost);
bestFitness = min(Cost);
best_pos = pop(index(1), :);

% Initialize Best Cost History
BestCosts = zeros(Max_FEs, 1);

%% Ivy Main Loop


it = 1;

% for FEs = 1:MaxFEs
while FEs <= Max_FEs

    % Get Best and Worst Cost Values
    % Costs = [pop.Cost];
    % % Best and worst costs
    % BestCost = min(Costs);
    % WorstCost = max(Costs);

    % Initialize new Population
    newpop = [];

    for i = 1:N

        ii=i+1;
        if i==N
            ii=1;
        end

        beta_1=1+(rand/2); % beta value in line 8 of "Algorithm 1 (Fig.2 in paper)"

        if  Cost(i)<beta_1*bestFitness
            %% Eq.(5)-(6)
            newsol =pop(i, :) +abs(randn(1, dim)) .*(pop(ii, :) - pop(i, :) )+ randn(1, dim) .*(GV(i, :));
        else
            %% Eq.(7)
            newsol =best_pos .*(rand+randn(1, dim).*(GV(i, :)));
        end
        
        %% Eq.(3)
        GV(i, :)=(GV(i, :)).*((rand^2)*(randn(1,dim)));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        newsol = max(newsol, lb);
        newsol = min(newsol, ub);
        %% Eq.(8)
        newGV(i, :)=(newsol /(ub-lb));

        %% Evaluate new population
        newCost(i) = fobj(newsol);
        FEs = FEs + 1;

        newpop = [newpop;
                  newsol];

    end

    % Merge Populations
    AllCost = [Cost; newCost];
    pop = [pop;
        newpop];

    % Sort Population
    [SortCost, SortOrder]=sort(AllCost);
    pop = pop(SortOrder, :);
    Cost = SortCost(1:N);
    % Competitive Exclusion (Delete Extra Members)
    pop = pop(1:N, :);
    Cost = Cost(1:N);

    % Store Best Solution Ever Found
    BestSol = pop(1);
    BestCosts(FEs) = Cost(1);
    Convergence_curve(it) =Cost(1);
    it = it + 1;
end
end