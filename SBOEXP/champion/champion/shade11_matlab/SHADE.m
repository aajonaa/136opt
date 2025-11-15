%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB/Octave source code of SHADE 1.1.
%% Please note that this source code is transferred from the C++ source code version.
%% About SHADE 1.1, please see following papers:
%%
%% * Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
%%
%% For this package, we downloaded JADE's source code from Dr. Q. Zhang's website (http://dces.essex.ac.uk/staff/qzhang) and modified it.
%%
%% Update
%% 9/Oct/2014: incompatibilities-fix between Octave and Matlab  (thanks to Dr. Elsayed)

%% this codes are slightly modified to cater for CEC 2005 benchmarks
%%%%%%%%%%%%%%%%%%%
function [GBEST, cg_curve]=GL_25(N,MaxFEs,lb,ub,dim,fobj)

format long;
format compact;

GBEST = zeros(1,dim);
fes = 0;
cg_curve=[];
children_fitness = zeros(N, 1);
t = 1;
pop_size = N;
fitness = zeros(N, 1);
problem_size = dim;
max_nfes = 10000 * problem_size;

rand('seed', sum(100 * clock));

lu = [lb .* ones(1, problem_size); ub .* ones(1, problem_size)];


    

%%  parameter settings for SHADE
p_best_rate = 0.1;
arc_rate = 2;
memory_size = problem_size;


%% Initialize the main population
popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
pop = popold; % the old population becomes the current population

%         fitness = feval(fhd,pop',func);
for i = 1 : N
    fitness(i) = fobj(pop(i, :));
    fes = fes + 1;
end
%         fitness = fitness';

nfes = 0;
bsf_fit_var = 1e+30;

%%%%%%%%%%%%%%%%%%%%%%%% for out
for i = 1 : pop_size
    nfes = nfes + 1;

    if fitness(i) < bsf_fit_var
        bsf_fit_var = fitness(i);
        bsf_solution = pop(i, :);
    end

    %% if mod(nfes, 1000) == 0
    %% bsf_error_var = bsf_fit_var  - optimum;
    %% if bsf_error_var < epsilon; bsf_error_var = 0; end;
    %%	fprintf(sprintf('%1.16e \n', bsf_error_var));
    %%    fprintf(sprintf('%d %1.16e \n', nfes, bsf_error_var));
    %% end

    %%      if nfes > max_nfes; exit(1); end
    % if nfes > max_nfes; break; end
end
%%%%%%%%%%%%%%%%%%%%%%%% for out

memory_sf = 0.5 .* ones(memory_size, 1);
memory_cr = 0.5 .* ones(memory_size, 1);
memory_pos = 1;

archive.NP = arc_rate * pop_size; % the maximum size of the archive
archive.pop = zeros(0, problem_size); % the solutions stored in te archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions

%% main loop
while fes <= MaxFEs

%     'SHADE'
%     t

    pop = popold; % the old population becomes the current population
    [~, sorted_index] = sort(fitness, 'ascend');

    mem_rand_index = ceil(memory_size * rand(pop_size, 1));
    mu_sf = memory_sf(mem_rand_index);
    mu_cr = memory_cr(mem_rand_index);

    %% for generating crossover rate
    cr = normrnd(mu_cr, 0.1);
    term_pos = find(mu_cr == -1);
    cr(term_pos) = 0;
    cr = min(cr, 1);
    cr = max(cr, 0);

    %% for generating scaling factor
    sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
    pos = find(sf <= 0);

    while ~ isempty(pos)
        sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
        pos = find(sf <= 0);
    end

    sf = min(sf, 1);

    r0 = [1 : pop_size];
    popAll = [pop; archive.pop];
    [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);

    pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
    randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
    randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
    pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions

    vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));


    vi = boundConstraint_SHADE(vi, pop, lu);

    mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
    rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
    jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
    ui = vi; ui(mask) = pop(mask);

%             children_fitness = feval(fhd, ui', func);
    for i = 1 : N
        children_fitness(i) = fobj(ui(i, :));
        fes = fes + 1;
    end
%             children_fitness = children_fitness';

    %%%%%%%%%%%%%%%%%%%%%%%% for out
    for i = 1 : pop_size
        nfes = nfes + 1;

        if children_fitness(i) < bsf_fit_var
            bsf_fit_var = children_fitness(i);
            bsf_solution = ui(i, :);
        end

        %% if mod(nfes, 1000) == 0
        %% bsf_error_var = bsf_fit_var  - optimum;
        %% if bsf_error_var < epsilon; bsf_error_var = 0; end;
        %%       fprintf(sprintf('%1.16e \n', bsf_error_var));
        %%       fprintf(sprintf('%d %1.16e \n', nfes, bsf_error_var));
        %%end

        %%	if nfes > max_nfes; exit(1); end
        % if nfes > max_nfes; break; end
    end
    %%%%%%%%%%%%%%%%%%%%%%%% for out

    dif = abs(fitness - children_fitness);


    %% I == 1: the parent is better; I == 2: the offspring is better
    I = (fitness > children_fitness);
    goodCR = cr(I == 1);
    goodF = sf(I == 1);
    dif_val = dif(I == 1);

    %      isempty(popold(I == 1, :))
    
    archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));

    [fitness, I] = min([fitness, children_fitness], [], 2);

    popold = pop;
    popold(I == 2, :) = ui(I == 2, :);

    num_success_params = numel(goodCR);

    if num_success_params > 0
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;

        %% for updating the memory of scaling factor
        memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);

        %% for updating the memory of crossover rate
        if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
            memory_cr(memory_pos)  = -1;
        else
            memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
        end

        memory_pos = memory_pos + 1;
        if memory_pos > memory_size;  memory_pos = 1; end
    end

    cg_curve(t) = bsf_fit_var;
    t = t + 1;
end
