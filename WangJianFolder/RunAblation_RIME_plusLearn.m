% Ablation on RIME_plusLearn components
clear; clc; close all;

SearchAgents_no = 30;
MaxFEsPerDim = 5000;
dim = 30;
Fold = 4;
FunctionSeq = {'F1','F2','F3','F4','F5'};
Function_name_all = {'F107','F108','F109','F110','F111'};

algorithmName = {
    'RIME_plusLearn', ...
    'RIME_plusLearn_noSF', ...
    'RIME_plusLearn_noD4', ...
    'RIME_plusLearn_smallElite', ...
    'RIME_plusLearn_bigTopK', ...
    'RIME_plusLearn_noStagnation' ...
};

MaxFEs = SearchAgents_no * MaxFEsPerDim;
foldNum = numel(algorithmName);
funcNum = numel(Function_name_all);
bestFitnessAll = inf * ones(foldNum, Fold, funcNum);

alg_fhds = cell(1, foldNum);
for a = 1:foldNum, alg_fhds{a} = str2func(algorithmName{a}); end

for func = 1:funcNum
    [lb,ub,~,fhd] = Get_Functions(Function_name_all{func}, dim);
    fprintf('---- RIME_plusLearn Ablation on %s ----\n', FunctionSeq{func});
    parfor fold = 1:Fold
        localBest = zeros(foldNum,1);
        for a = 1:foldNum
            [~, curve] = alg_fhds{a}(SearchAgents_no, MaxFEs, lb, ub, dim, fhd);
            localBest(a,1) = curve(end);
        end
        bestFitnessAll(:, fold, func) = localBest;
    end
end

AvgFitness = squeeze(mean(bestFitnessAll, 2));
AvgFitnessRank = zeros(size(AvgFitness));
for col = 1:funcNum, [~,~,r] = unique(AvgFitness(:,col)); AvgFitnessRank(:,col) = r; end
AvgFitnessRankSum = sum(AvgFitnessRank, 2);

disp('AvgFitnessRank per function:');
for a = 1:foldNum
    fprintf('%-28s : %s\n', algorithmName{a}, num2str(AvgFitnessRank(a,:)));
end
disp('AvgFitnessRankSum and AvgRank:');
for a = 1:foldNum
    fprintf('%-28s : sum=%d  avg=%.3f\n', algorithmName{a}, AvgFitnessRankSum(a), AvgFitnessRankSum(a)/funcNum);
end

