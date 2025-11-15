%% Use to benckmark the CEC 2017 functions
% function FitnessAll = ParameterAnalysis(N, Max_FEs, lb, ub, dim)

% general parameters
N = 30;
Max_FEs = 300000;
lb = -100;
ub = 100;
dim = 30;

% custom parameters
funcNum = 30; % CEC 2017
paraNum = 8; %k from 1/2 -- 1/9

%% Test 30 function and get the best fitness of 30 functions
% fitness matrix, every row is a result of corresponding parameter k
% every row have 30 columns which represent the result of corresponding
% paramter k
FitnessAll = ones(paraNum, funcNum) * inf; 
bestFitness = [];

for temp = 2:9
    disp(temp - 1);
    k = 1/temp
    
    % test 30 functions of CEC 17 for every parameter k
    % Simulate flod = 4
    for i = 1:funcNum
        fobj = @(x) cec17_func(x', i);
        bestFitness(i) = AOPara(N, Max_FEs, lb, ub, dim, fobj, k);
    end
    %% k is the parameter in the AO algorithm
    FitnessAll(temp-1, :) = bestFitness;
    
end

ranked_FitnessAll = ones(size(FitnessAll)) * inf;

for col = 1:size(FitnessAll, 2)
    [~, ~, rank] = unique(FitnessAll(:, col));  % Get the rank of each element in the column
    ranked_FitnessAll(:, col) = rank;  % Store the rank in the corresponding column
end

% Sum the every row value
row_sum_ranked_FitnessAll = sum(ranked_FitnessAll, 2);

[~, bestRow] = min(row_sum_ranked_FitnessAll);

% Display the results
disp('Original matrix:');
disp(FitnessAll);
disp('Ranked matrix:');
disp(ranked_FitnessAll);
disp('Sum of ranks for each row:');
disp(row_sum_ranked_FitnessAll);
disp(['The best row is: Row ', num2str(bestRow)]);
