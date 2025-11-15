%**************************************************************************************************
%Reference:  N. Hansen and A. Ostermeier, ¡°Completely derandomized self-adaptation
%                     in evolution strategies,¡± Evolut. Comput., vol. 9, no. 2, pp. 159¨C195, 2001.
%
% Note: We obtained the MATLAB source code from the authors, and did some
%           minor revisions in order to solve the 25 benchmark test functions,
%           however, the main body was not changed.
%**************************************************************************************************

function [GBEST,cg_curve]=CMAES(N,MaxFEs,lb,ub,dim,fobj)

tic;

format long;
format compact;

bestval = inf;
fes = 0;
GBEST = zeros(1,dim);
cg_curve=[];
t = 1;
lu = [lb * ones(1, dim); ub * ones(1, dim)];

% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.

% Define the dimension of the problem
n = dim;

rand('seed', sum(100 * clock));

% Initialize the main population

N = n;
maxeval = 300000;  % stop criteria
xmeanw = (lu(1, :) + rand(1, N) .* (lu(2, :) - lu(1, :)))'; % object parameter start point
sigma = 0.25; minsigma = 1e-15; maxsigma = max(lu(2, :)-lu(1, :)) / sqrt(n); % initial step size, minimal step size
flginiphase = 1; % Initial phase

% Parameter setting: selection,
lambda = 4 + floor(3 * log(N)); mu = floor(lambda/2);

arweights = log((lambda + 1)/2) - log(1:mu)'; % muXone array for weighted recomb.
% lambda = 10; mu = 2; arweights = ones(mu, 1); % uncomment for (2_I, 10)-ES
% parameter setting: adaptation
cc = 4/(N + 4); ccov = 2/(N + 2^0.5)^2;
cs = 4/(N + 4); damp = (1 - min(0.7, N * lambda/maxeval)) / cs + 1;

% Initialize dynamic strategy parameters and constants
B = eye(N); D = eye(N); BD = B * D; C = BD * transpose(BD);
pc = zeros(N, 1); ps = zeros(N, 1);
cw = sum(arweights) / norm(arweights); chiN = N^0.5 * (1 - 1/(4 * N) + 1/(21 * N^2));

% Generation loop
%     disp(['  (' num2str(mu) ', ' num2str(lambda) ')-CMA-ES (w = [' num2str(arweights', '%5.2f') '])' ]);
counteval = 0; flgstop = 0;

% Boundary
lb = (ones(lambda, 1) * lu(1, :))';
ub = (ones(lambda, 1) * lu(2, :))';

while fes <= MaxFEs % if flgstop break;

    
    % Generate and evaluate lambda offspring
    arz = randn(N, lambda);
    arx = xmeanw * ones(1, lambda) + sigma * (BD * arz);
    x_ = xmeanw * ones(1, lambda);

    % You may handle constraints here and now.
    % You may either resample columns of arz and/or multiply
    % them with a factor < 1 (the latter will result in a
    % decreased overall step size) and recalculate arx
    % accordingly. Do not change arx or arz in any other
    % way. You may also use an altered arx only for the
    % evaluation of the fitness function, if you leave
    % arx and arz unchanged for the algorithm.

    % Handle the elements of the variable which violate the boundary
    I = find(arx > ub);
    arx(I) = 2 * ub(I) - arx(I);
    aa = find(arx(I) < lb(I));
    arx(I(aa)) = lb(I(aa));
    I = find(arx < lb);
    arx(I) = 2 * lb(I) - arx(I);
    aa = find(arx(I) > ub(I));
    arx(I(aa)) = ub(I(aa));

    U = arx';
    arfitness = zeros(lambda, 1);
    for i = 1 : lambda
        arfitness(i) = fobj(U(i, :));
        fes = fes + 1;
    end
    counteval = counteval + lambda;
    % Sort by fitness and compute weighted mean in xmeanw
    [arfitness, arindex] = sort(arfitness); % minimization
    xold = xmeanw; % for speed up of Eq. (14)
    xmeanw = arx(:, arindex(1:mu)) * arweights/sum(arweights);
    zmeanw = arz(:, arindex(1:mu)) * arweights/sum(arweights);

    % Adapt covariance matrix
    pc = (1-cc) * pc + (sqrt(cc * (2-cc)) * cw/sigma) * (xmeanw-xold); % Eq. (14)
    if ~flginiphase % do not adapt in the initial phase
        C = (1-ccov) * C + ccov * pc * transpose(pc);           % Eq. (15)
    end
    % adapt sigma
    ps = (1-cs) * ps + (sqrt(cs * (2-cs)) * cw) * (B * zmeanw);      % Eq. (16)
    sigma = sigma * exp((norm(ps)-chiN)/chiN/damp);        % Eq. (17)

    % Update B and D from C
    if mod(counteval/lambda, 1/ccov/N/5) < 1
        C = triu(C) + transpose(triu(C, 1)); % enforce symmetry
        [B, D] = eig(C);
        % limit condition of C to 1e14 + 1
        if max(diag(D)) > 1e14 * min(diag(D))
            tmp = max(diag(D))/1e14 - min(diag(D));
            C = C + tmp * eye(N); D = D + tmp * eye(N);
        end
        D = diag(sqrt(diag(D))); % D contains standard deviations now
        BD = B * D; % for speed up only
    end % if mod

    % Adjust minimal step size
    if sigma * min(diag(D)) < minsigma ...
            | arfitness(1) == arfitness(min(mu + 1, lambda)) ...
            | xmeanw == xmeanw ...
            + 0.2 * sigma * BD(:, 1 + floor(mod(counteval/lambda, N)))
        sigma = 1.4 * sigma;

        % flgstop = 1;
    end
    if sigma > maxsigma
        sigma = maxsigma;
    end

    % Test for end of initial phase
    if flginiphase && counteval/lambda > 2/cs
        if (norm(ps)-chiN)/chiN < 0.05 % step size is not much too small
            flginiphase = 0;
        end
    end
    
    if ( arfitness(1) < bestval )
        bestval = arfitness(1);
    end
    cg_curve(t) = bestval;
    t = t + 1;
end % while, end generation loop

%     disp([num2str(counteval) ': ' num2str(arfitness(1))]);
%     if exist('sfile', 'var')
%       disp(['Results saved in ' sfile]);
%     end

toc;
end