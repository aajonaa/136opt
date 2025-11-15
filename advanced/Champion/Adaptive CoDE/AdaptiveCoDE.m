%**************************************************************************************************
%Author: Yong Wang
%Last Edited: July 1, 2010
%Email: ywang@csu.edu.cn
%Reference: Differential Evolution with Composite Trial Vector Generation Strategies
%                                                       and Control Parameters
%                            IEEE Transactions on Evolutionary Computation, Accepted
%**************************************************************************************************

function [GBEST, cg_curve]=GL_25(N,MaxFEs,lb,ub,dim,fobj)
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

% Learning period
LEP = 50;

rand('seed', sum(100 * clock));

% Initialize the main population
p = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));

% Evaluate the objective function values
for i = 1 : N
    fit(i) = fobj(p(i, :));
    fes = fes + 1;
end

for j = 1 : 3

    % Success memory for each trial vector generation strategy
    sucMemo{ j } = zeros(LEP, 3);
    % Failure memory for each trial vector generation strategy
    failMemo{ j } = zeros(LEP, 3);

end

% The number of function evaluations (FES)
FES = popsize;

gen = 1;

% while FES < n * 10000
while fes <= MaxFEs

    'AdaptiveCoDE'
    t

    if gen < LEP

        % The control parameter settings have the same probability
        % to be selected for each trial vector generation strategy.
        paraProb = ones(3, 3) * 1/3;

    else

        for k = 1 : 3

            % Compute the success rate of each control parameter
            % setting for each trial vector generation strategy
            paraSR = sum(sucMemo{ k }) ./ (sum(sucMemo{ k }) + sum(failMemo{ k })) + 0.01;
            % Normalized the success rates
            paraProb(k, :) = paraSR ./ sum(paraSR);

        end

    end

    pTemp = p;
    fitTemp = fit;

    % uSet: the set of trial vectors
    uSet = zeros(3 * popsize, n);

    % numP{ i }: at each generation, record which control parameter
    % setting is used to produce the trial vector for each target
    % vector, please note that three trial vectors are produced for
    % each target vector
    numP = cell(popsize, 1);
    for i = 1 : popsize
        numP{ i } = zeros(3);
    end
    % numS: at each generation, record which control parameter
    % setting is used to produce the trial vector entering the next
    % population successfully
    numS = zeros(3);

    for i = 1 : popsize

        % Generate the trail vectors
        [u, tempPara] = reproduce(p, lu, i, popsize, n, paraProb);

        uSet(i * 3-2 : 3 * i, :) = u;

        FES = FES + 3;

        for k = 1 : 3

            % Judge which control parameter setting is used to
            % produce the trial vector for each target vector,
            % please note that three trial vectors are generated
            numP{ i }(k, 1) = numP{ i }(k, 1) + length(find(tempPara(k, 1) == 1));
            numP{ i }(k, 2) = numP{ i }(k, 2) + length(find(tempPara(k, 1) == 2));
            numP{ i }(k, 3) = numP{ i }(k, 3) + length(find(tempPara(k, 1) == 3));

        end

    end

    % Evaluate the trial vectors
    fitSet = zeros(3*N, 1);
    for i = 1 : 3 * popsize
        fitSet(i) = fobj(uSet(i, :));
        fes = fes + 1;
    end

    for i = 1 : popsize

        % Find the best trial vector of the three trial vectors
        % generated for each target vector
        % Minindex: denotes which trial vector generation strategy
        % is used to generate the best trial vector
        [~, minIndex] = min(fitSet(3 * i - 2 : 3 * i, :));
        bestInd = uSet(3 * (i - 1) + minIndex, :);
        bestIndFit = fitSet(3 * (i - 1) + minIndex, :);

        % Selection
        if fit(i) >= bestIndFit

            pTemp(i, :) = bestInd;
            fitTemp(i, :) = bestIndFit;

            % Judge which control parameter setting is used to
            % generate the best trial entering the next population
            temp = find(numP{ i }(minIndex, :) == 1);
            numS(minIndex, temp) = numS(minIndex, temp) + 1;

        end

    end

    p = pTemp;
    fit = fitTemp;

    % Update the success memory
    for k = 1 : 3

        sucMemo{ k }(1, :) = [];
        sucMemo{ k }(LEP, :) = numS(k, :);

    end

    % Record the total number of each control parameter setting
    % used at each generation for each trial vector generation
    % strategy
    totalNum = zeros(3);
    for i = 1 : popsize

        totalNum = totalNum + numP{ i };

    end

    % Update the failure memory
    for k = 1 : 3

        failMemo{ k }(1, :) = [];
        failMemo{ k }(LEP, :) = totalNum(k, :) - numS(k, :);

    end

    gen = gen + 1;



    cg_curve(t) = min(fit);
    t = t + 1;
end

toc;
