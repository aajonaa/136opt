%**************************************************************************************************
%Reference:  1) J. Brest, S. Greiner, B. Boskovic, M. Mernik, and V. Zumer, ¡°Self-adapting
%                     control parameters in differential evolution: A comparative study on numerical
%                     benchmark problems,¡± IEEE Trans. Evolut. Comput., vol. 10, no. 6,
%                     pp. 646¨C657, Dec. 2006.
%                     2) J. Zhang and A. C. Sanderson, ¡°JADE: adaptive differential evolution
%                     with optional external archive,¡± IEEE Trans. Evolut. Comput., vol. 13,
%                     no. 5, pp. 945-958, 2009.
%
% Note: We obtained the MATLAB source code from Dr. J. Zhang, and did some
%           minor revisions in order to solve the 25 benchmark test functions,
%           however, the main body was not changed.

% This modified Matlab codes are downloaded from http://ist.csu.edu.cn/YongWang.htm
%**************************************************************************************************

function [GBEST, convergence]=JDE(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.
if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
% Define the dimension of the problem
n = dim;
lu=[lb;ub];

GBEST = zeros(1,dim);
fes = 0;
popsize = SearchAgents_no;

tau1 = 0.1;
tau2 = 0.1;

F = 0.5 * ones(popsize, 1);
CR = 0.9 * ones(popsize, 1);

% convergence = [];  % record the best results

%Main body which was provided by Dr. J. Zhang

rand('seed', sum(100 * clock));

% Initialize the main population
popold = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :)-lu(1, :), popsize, 1));

valParents=zeros(popsize,1);
for i=1:popsize
    valParents(i) = fobj(popold(i,:));
    fes = fes + 1;
end

% FES = 0;
l=1;
convergence = [];
% while FES < n * 10000
while fes<MaxFEs    
    pop = popold;      % the old population becomes the current population
    
    Fold = F;
    CRold = CR;
    
    IF  = rand(popsize, 1) < tau1;
    ICR = rand(popsize, 1) < tau2;
    
    F(IF)  = 0.1 + 0.9 * rand(sum(IF), 1);
    CR(ICR) = 0.0 + 1.0 * rand(sum(ICR), 1);
    
    r0 = [1:popsize];
    [r1, r2, r3] = gnR1R2R3(popsize, r0);
    
    %  == == == == == = Mutation == == == == == == == == =
    vi  = pop(r1, :) + F(:, ones(1, n)) .* (pop(r2, :) - pop(r3, :));
    
    vi = boundConstraint_JDE(vi, lu);
    
    % == == == == =  Crossover == == == == =
    
    mask = rand(popsize, n) > CR(:, ones(1, n));     % mask is used to indicate which elements of ui comes from the parent
    rows = (1:popsize)'; cols = floor(rand(popsize, 1) * n) + 1;  % choose one position where the element of ui doesn't come from the parent
    jrand = sub2ind([popsize n], rows, cols); mask(jrand) = false;
    ui = vi;  ui(mask) = pop(mask);
    
    valOffspring=zeros(popsize,1);
    for i=1:popsize
        valOffspring(i) = fobj(ui(i,:));
        fes = fes + 1;
    end
    
%     FES = FES + popsize;
    
    %  == == == == == = Selection == == == == == =
    % I == 1: the parent is better; I == 2: the offspring is better
    [valParents, I] = min([valParents, valOffspring], [], 2);
    popold = pop;
    
    popold(I == 2, :) = ui(I == 2, :);
    
    F(I == 1) = Fold(I == 1);
    CR(I == 1) = CRold(I == 1);
    
%     convergence = [convergence min(valParents)];
      convergence(l)=min(valParents);
      l=l+1;    
end

% bestScore=convergence(end);

end

