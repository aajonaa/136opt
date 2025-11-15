% Ablation runner for RIME, GO, and RIME_GO variants
clear; clc; close all;

% User knobs
SearchAgents_no = 30;
MaxFEsPerDim = 5000;    % baseline from your main script
dim = 30;               % CEC17 default
Fold = 4;               % repetitions per function
FunctionSeq = {'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10'};
Function_name_all = {'F107','F108','F109','F110','F111','F112','F113','F114','F115','F116'};

% Choose algorithms for ablation
algorithmName = {
    'RIME', ...
    'GO', ...
    'RIME_GO', ...
    'RIME_GO_noGO', ...        % isolate RIME part
    'RIME_GO_onlyGO', ...      % isolate GO part
    'RIME_GO_noPrune', ...     % prune effect
    'RIME_GO_noAdaptiveN', ... % adaptive N effect
    'RIME_GO_onlyLearn', ...   % GO learning half
    'RIME_GO_onlyReflect' ...  % GO reflection half
};

MaxFEs = SearchAgents_no * MaxFEsPerDim; % consistent with main

foldNum = numel(algorithmName);
funcNum = numel(Function_name_all);
bestFitnessAll = inf * ones(foldNum, Fold, funcNum);

% Prepare handles
alg_fhds = cell(1, foldNum);
for a = 1:foldNum
    alg_fhds{a} = str2func(algorithmName{a});
end

for func = 1:funcNum
    Function_name = Function_name_all{func};
    [lb,ub,~,fhd] = Get_Functions(Function_name, dim);
    fprintf('---- Ablation on %s (%s) ----\n', FunctionSeq{func}, Function_name);

    parfor fold = 1:Fold
        localBest = zeros(foldNum,1);
        for a = 1:foldNum
            [~, curve] = alg_fhds{a}(SearchAgents_no, MaxFEs, lb, ub, dim, fhd);
            localBest(a,1) = curve(end);
        end
        bestFitnessAll(:, fold, func) = localBest;
    end
end

% Aggregate
AvgFitness = squeeze(mean(bestFitnessAll, 2));
StdFitness = squeeze(std(bestFitnessAll, 0, 2));

% Ranks per function
AvgFitnessRank = zeros(size(AvgFitness));
for col = 1:funcNum
    [~,~,r] = unique(AvgFitness(:,col));
    AvgFitnessRank(:,col) = r;
end
AvgFitnessRankSum = sum(AvgFitnessRank, 2);
[~,~,AlgoRank] = unique(AvgFitnessRankSum);

% Print summary
disp('AvgFitnessRank per function:');
for a = 1:foldNum
    fprintf('%-20s : %s\n', algorithmName{a}, num2str(AvgFitnessRank(a,:)));
end

disp('AvgFitnessRankSum and AvgRank:');
for a = 1:foldNum
    avgRank = AvgFitnessRankSum(a) / funcNum;
    fprintf('%-20s : sum=%d  avg=%.3f\n', algorithmName{a}, AvgFitnessRankSum(a), avgRank);
end

% Optional: write XLS (disabled by default)
% xlswrite logic can be added similarly to ImproveTest_v2 if needed.

