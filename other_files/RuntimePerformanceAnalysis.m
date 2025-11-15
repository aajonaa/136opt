clear all;
close all;
clc;

% Define algorithm names
% algorithmName = {'SBO', 'TLBO', 'PO', 'HHO', 'RIME', 'ALCPSO', 'CLPSO', 'CGPSO', 'MSPSO', 'CMAES', 'DECLS', 'LSHADE_cnEpSin', 'SCADE', 'SHADE'};
% algorithmName = {'DECLS', 'SBO'};
algorithmName = {'LSHADE_cnEpSin'};

% Initialize a table to store execution times
% dimensions = [10];
dimensions = [10, 30, 50, 100];
executionTimes = zeros(length(dimensions), length(algorithmName)); % Rows: dimensions, Columns: algorithms

% Loop through dimensions
for dimIdx = 1:length(dimensions)
    dim = dimensions(dimIdx);
    fprintf('Running for dimension: %d\n', dim);
    
    % Loop through algorithms
    for algo = 1:length(algorithmName)
        algo_fhd = str2func(algorithmName{algo}) % Get function handle for the algorithm
        
        tic; % Start timer
        for funcNum = 1:30
            if funcNum == 2 % Skip F108
                continue;
            end
            
            fhd = @(x) cec17_func(x', funcNum); % Define the objective function
            parfor cflod = 1:4 % Run 30 times for each function
                [~, ~] = algo_fhd(30, 300000, -100, 100, dim, fhd); % Run the algorithm
            end
        end
        elapsedTime = toc % Stop timer and get elapsed time
        executionTimes(dimIdx, algo) = elapsedTime; % Store execution time
    end
end

% Create a table for the results
resultsTable = array2table(executionTimes, 'VariableNames', algorithmName, 'RowNames', string(dimensions));
resultsTable.Properties.DimensionNames{1} = 'Dimension'; % Name the first dimension
resultsTable.Properties.DimensionNames{2} = 'Algorithm'; % Name the second dimension

% Display the table
disp(resultsTable);

% Write the table to an Excel file
% outputFileName = 'AlgorithmExecutionTimes.xlsx';
outputFileName = 'LSHADE_cnEpSinExecutionTimes.xlsx';
writetable(resultsTable, outputFileName, 'WriteRowNames', true);

fprintf('Execution times have been written to %s\n', outputFileName);