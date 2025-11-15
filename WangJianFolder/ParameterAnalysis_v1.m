%% Use to benckmark the CEC 2017 functions With single function but different parameters

%% Global parameters
N = 30;
Max_FEs = 300000;
lb = -100;
ub = 100;
dim = 30;

%% custom parameters
funcNum = 30; % CEC 2017

% %% Parameter serve as different function
% paraNum = 3;
% para = {'sin', 'cos', 'tan'};
% para_handle = cell(1, paraNum);
% for i = 1:paraNum
%     para_handle{i} = str2func(para{i});
% end
%% Parameter serve as constant
paraNum = 9;

%% Test 30 function and get the best fitness of 30 functions
FitnessAll = ones(paraNum, funcNum) * inf; 
bestFitness = ones(1, funcNum) * inf;

%% Row: Paramter Col: Fitness
for para = 1:paraNum
%     %% Here is where parameter need to send to algorithm for function
%     para_send = para_handle{para} %#ok<NOPTS>
    %% Here is where parameter need to send to algorithm for constant
    para_send = 1/10 * para %#ok<NOPTS>
    bestFitnessLocal = inf * ones(1, paraNum);
    % Bench mark every single function here
    for func = 1:funcNum
        if func == 2
            continue;
        end
        Function_name=['F',num2str(func)];
        disp(['----------------',Function_name,'--------------------']);
        fobj = @(x) cec17_func(x', func);
        [best_pos, ~] = AOPara2(N, Max_FEs, lb, ub, dim, fobj, para_send);
        bestFitnessLocal(func) = fobj(best_pos);
    end
    FitnessAll(para, :) = bestFitnessLocal;
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
disp(['The best para is: Row ', num2str(bestRow)]);

