% Test script to verify all dependencies are present and the folder is self-contained
% Run this script from the ExplorExploi folder to verify everything is working

fprintf('=== Testing ExplorExploi Folder Dependencies ===\n\n');

% Test 1: Check if main script exists
fprintf('Test 1: Checking main script... ');
if exist('exp_CCMWOA1.m', 'file')
    fprintf('PASS\n');
else
    fprintf('FAIL - exp_CCMWOA1.m not found\n');
end

% Test 2: Check algorithm files
fprintf('Test 2: Checking algorithm files... ');
if exist('LRRIME/LRRIME.m', 'file') && exist('RMRIME/RIME.m', 'file')
    fprintf('PASS\n');
else
    fprintf('FAIL - Algorithm files not found\n');
end

% Test 3: Check utility functions
fprintf('Test 3: Checking utility functions... ');
utils = {'subFunction/Get_Functions.m', 'subFunction/MySampling.m', ...
         'subFunction/initialization.m', 'Divergence_calculation.m'};
all_exist = true;
for i = 1:length(utils)
    if ~exist(utils{i}, 'file')
        all_exist = false;
        fprintf('\n  Missing: %s', utils{i});
    end
end
if all_exist
    fprintf('PASS\n');
else
    fprintf('\nFAIL\n');
end

% Test 4: Check statistical functions
fprintf('Test 4: Checking statistical functions... ');
stats = {'statics/Orderhao.m', 'statics/pValueToExcelhao.m', ...
         'statics/FridTest3.m', 'statics/FridTest4.m'};
all_exist = true;
for i = 1:length(stats)
    if ~exist(stats{i}, 'file')
        all_exist = false;
        fprintf('\n  Missing: %s', stats{i});
    end
end
if all_exist
    fprintf('PASS\n');
else
    fprintf('\nFAIL\n');
end

% Test 5: Check CEC benchmark functions
fprintf('Test 5: Checking CEC benchmark functions... ');
cec_funcs = {'cec05_func.m', 'cec13_func.mexw64', 'cec14_func.mexw64', ...
             'cec17_func.mexw64', 'cec19_func.mexw64'};
all_exist = true;
for i = 1:length(cec_funcs)
    if ~exist(cec_funcs{i}, 'file')
        all_exist = false;
        fprintf('\n  Missing: %s', cec_funcs{i});
    end
end
if all_exist
    fprintf('PASS\n');
else
    fprintf('\nFAIL\n');
end

% Test 6: Check data folders
fprintf('Test 6: Checking data folders... ');
if exist('input_data_cec2017', 'dir') && exist('cec05', 'dir')
    fprintf('PASS\n');
else
    fprintf('FAIL - Data folders not found\n');
end

% Test 7: Check output folder
fprintf('Test 7: Checking output folder... ');
if exist('exp_result', 'dir')
    fprintf('PASS\n');
else
    fprintf('FAIL - exp_result folder not found\n');
end

% Test 8: Test basic function calls
fprintf('Test 8: Testing basic function calls... ');
try
    % Test initialization
    test_pop = initialization(5, 10, 100, -100);
    
    % Test Get_Functions
    [lb, ub, dim, fobj] = Get_Functions('F1', 10);
    
    % Test MySampling
    test_data = rand(1, 100);
    sampled = MySampling(test_data, 10);
    
    % Test Divergence_calculation
    div = Divergence_calculation(test_pop, 10, 5);
    
    fprintf('PASS\n');
catch ME
    fprintf('FAIL\n');
    fprintf('  Error: %s\n', ME.message);
end

% Test 9: Verify CEC2017 can be called
fprintf('Test 9: Testing CEC2017 function... ');
try
    test_x = zeros(30, 1);
    result = cec17_func(test_x, 1);
    fprintf('PASS\n');
catch ME
    fprintf('FAIL\n');
    fprintf('  Error: %s\n', ME.message);
end

fprintf('\n=== All Tests Complete ===\n');
fprintf('If all tests passed, the folder is self-contained and ready to use.\n');
fprintf('You can now run exp_CCMWOA1.m without errors.\n');

