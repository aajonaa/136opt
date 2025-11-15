%**************************************************************************************************
%Author: Yong Wang
%Last Edited: July 1, 2010
%Email: ywang@csu.edu.cn
%Reference: Differential Evolution with Composite Trial Vector Generation Strategies
%                                                    and Control Parameters
%                           IEEE Transactions on Evolutionary Computation, Accepted
%**************************************************************************************************

function [GBEST, cg_curve]=CoDE(N,MaxFEs,lb,ub,dim,fobj)
tic;

format long;
format compact;

fes = 0;
GBEST = zeros(1,dim);
cg_curve=[];
t = 1;
fit = zeros(N, 1);
lu = [lb * ones(1, dim); ub * ones(1, dim)];

% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.

% Define the dimension of the problem
n = dim;

% Main body
popsize = N;

rand('seed', sum(100*clock));

% Initialize the main population
p = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));

% Evaluate the objective function values
for i = 1 : N
    fit(i) = fobj(p(i, :));
    fes = fes + 1;
end

% Record the number of function evaluations (FES)
FES = popsize;
% while FES < n * 10000
while fes <= MaxFEs
    
%     'CoDE'
%     t

    pTemp = p;
    fitTemp = fit;

    % uSet: the set of trial vectors
    uSet = zeros(3 * popsize, n);

    for i = 1 : popsize

        % the three control parameter settings
        F    = [1.0 1.0 0.8];
        CR = [0.1 0.9 0.2];

        % Uniformly and randomly select one of the control
        % parameter settings for each trial vector generation strategy
        paraIndex = floor(rand(1, 3) * length(F)) + 1;

        % Generate the trail vectors
        u = generator(p, lu, i, F, CR, popsize, n, paraIndex);

        uSet(i * 3 - 2 : 3 * i, :) = u;

        FES = FES + 3;

    end

    % Evaluate the trial vectors
    fitSet = zeros(3*N, 1);
    for i = 1 : 3 * popsize
        fitSet(i) = fobj(uSet(i, :));
        fes = fes + 1;
    end

    for i = 1 : popsize

        % Choose the best trial vector from the three trial vectors
        [~, minID] = min(fitSet(3 * i - 2 : 3 * i, :));
        bestInd = uSet(3 * (i - 1) + minID, :);
        bestIndFit = fitSet(3 * (i - 1) + minID, :);

        % Choose the better one between the trial vector and the
        % target vector
        if fit(i) >= bestIndFit

            pTemp(i, :) = bestInd;
            fitTemp(i, :) = bestIndFit;

        end

    end

    p = pTemp;
    fit = fitTemp;

    cg_curve(t) = min(fit);
    t = t + 1;
end
        
toc;
