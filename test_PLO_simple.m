%% Simple test to verify PLO algorithms are working correctly
clear; close all; clc;

% Add all paths
addpath(genpath(pwd));

% Test on a simple function first (Sphere function)
dim = 30;
lb = -100;
ub = 100;
MaxFEs = 30000; % Reduced for quick testing
N = 30;

% Simple sphere function
fobj = @(x) sum(x.^2);

fprintf('Testing PLO algorithms on Sphere function...\n');

% Test original PLO
fprintf('Testing PLO...\n');
try
    [best_pos1, conv_curve1] = PLO(N, MaxFEs, lb, ub, dim, fobj);
    fprintf('PLO final fitness: %e\n', conv_curve1(end));
catch ME
    fprintf('PLO failed: %s\n', ME.message);
end

% Test PLO_EPM_Strategy1
fprintf('Testing PLO_EPM_Strategy1...\n');
try
    [best_pos2, conv_curve2] = PLO_EPM_Strategy1(N, MaxFEs, lb, ub, dim, fobj);
    fprintf('PLO_EPM_Strategy1 final fitness: %e\n', conv_curve2(end));
catch ME
    fprintf('PLO_EPM_Strategy1 failed: %s\n', ME.message);
end

% Test PLO_EPM_Strategy2
fprintf('Testing PLO_EPM_Strategy2...\n');
try
    [best_pos3, conv_curve3] = PLO_EPM_Strategy2(N, MaxFEs, lb, ub, dim, fobj);
    fprintf('PLO_EPM_Strategy2 final fitness: %e\n', conv_curve3(end));
catch ME
    fprintf('PLO_EPM_Strategy2 failed: %s\n', ME.message);
end

% Plot convergence curves
figure;
if exist('conv_curve1', 'var')
    semilogy(conv_curve1, 'r-', 'LineWidth', 2); hold on;
end
if exist('conv_curve2', 'var')
    semilogy(conv_curve2, 'g-', 'LineWidth', 2); hold on;
end
if exist('conv_curve3', 'var')
    semilogy(conv_curve3, 'b-', 'LineWidth', 2); hold on;
end
xlabel('Iteration');
ylabel('Best Fitness');
title('PLO Variants Convergence on Sphere Function');
legend('PLO', 'PLO\_EPM\_Strategy1', 'PLO\_EPM\_Strategy2');
grid on;

fprintf('Test completed. Check the convergence plot.\n');
