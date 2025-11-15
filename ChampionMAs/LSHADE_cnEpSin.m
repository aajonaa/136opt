%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB/Octave source code of LSHADE_cnEpSin which is a new version of LSHADE-EpSin.
%% Please see the following papers:
%% 1. LSHADE_cnEpSin:
%%     Noor H. Awad, Mostafa Z. Ali, Ponnuthurai N. Suganthan, Ensemble Sinusoidal Differential Covariance Matrix Adaptation with Euclidean Neighborhood  for Solving CEC2017 Benchmark Problems, in Proc. IEEE Congr. Evol. Comput. CEC 2017, June, Donostia - San Sebastián, Spain

%% 2. LSHADE-EpSin:
%%    Noor H. Awad, Mostafa Z. Ali, Ponnuthurai N. Suganthan and Robert G. Reynolds: An Ensemble Sinusoidal Parameter Adaptation incorporated with L-SHADE for Solving CEC2014 Benchmark Problems, in Proc. IEEE Congr. Evol. Comput. CEC 2016, Canada, July, 2016

%% About L-SHADE, please see following papers:
%% Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
%%  J. Zhang, A.C. Sanderson: JADE: Adaptive differential evolution with optional external archive,?IEEE Trans Evol Comput, vol. 13, no. 5, pp. 945?58, 2009
function [best_pos,Convergence_curve] = LSHADE_cnEpSi(N,MaxFEs,lb,ub,dim,fobj)
tic
Convergence_curve = [];
format long;
format compact;
freq_inti = 0.5;
rand('seed', sum(100 * clock));
lu = [lb * ones(1, dim); ub * ones(1, dim)];
pb = 0.4;
ps = .5;
S.Ndim = dim;
S.Lband = ones(1, S.Ndim)*(-100);
S.Uband = ones(1, S.Ndim)*(100);
run_funcvals = [];     
%%  parameter settings for L-SHADE
p_best_rate = 0.11;    %0.11
arc_rate = 1.4;
memory_size = 5;
SEL = round(ps*N);
max_pop_size = N;
min_pop_size = 4.0;
FEs = 0;
%% Initialize the main population
popold = repmat(lu(1, :), N, 1) + rand(N, dim) .* (repmat(lu(2, :) - lu(1, :), N, 1));
X = popold; % the old population becomes the current population
for i=1:size(X,1)
    fitness(i)=fobj(X(i,:));
    FEs=FEs+1;
end
fitness = fitness';

[minFit, fitIndex] = sort(fitness);
bestFitness = minFit;
best_pos = X(fitIndex, :);

memory_sf = 0.5 .* ones(memory_size, 1);
memory_cr = 0.5 .* ones(memory_size, 1);

memory_freq = freq_inti*ones(memory_size, 1);
memory_pos = 1;

gg=0;  %%% generation counter used For Sin
goodF1all = [];
goodF2all =[];
badF1all = [];
badF2all = [];
goodF1 = [];
goodF2 = [];
l=1;
while FEs<MaxFEs
    gg=gg+1;
    X = popold; % the old population becomes the current population
    [~, sorted_index] = sort(fitness, 'ascend');
    
    mem_rand_index = ceil(memory_size * rand(N, 1));
    mu_sf = memory_sf(mem_rand_index);
    mu_cr = memory_cr(mem_rand_index);
    mu_freq = memory_freq(mem_rand_index);
    
    %% for generating crossover rate
    cr = normrnd(mu_cr, 0.1);
    term_pos = find(mu_cr == -1);
    cr(term_pos) = 0;
    cr = min(cr, 1);
    cr = max(cr, 0);
    
    %% for generating scaling factor
    sf = mu_sf + 0.1 * tan(pi * (rand(N, 1) - 0.5));
    pos = find(sf <= 0);
    while ~ isempty(pos)
        sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
        pos = find(sf <= 0);
    end
    freq = mu_freq + 0.1 * tan(pi*(rand(N, 1) - 0.5));
    pos_f = find(freq <=0);
    while ~ isempty(pos_f)
        freq(pos_f) = mu_freq(pos_f) + 0.1 * tan(pi * (rand(length(pos_f), 1) - 0.5));
        pos_f = find(freq <= 0);
    end
    sf = min(sf, 1);
    freq = min(freq, 1);
    LP = 20;
    flag1 = false;
    flag2 = false;
    if(FEs <= MaxFEs/2)
        flag1 = false;
        flag2 = false;
        if (gg <= LP)
            p1 = 0.5;
            p2 = 0.5;
            c=rand;
            if(c < p1)
                sf = 0.5.*( sin(2.*pi.*freq_inti.*gg+pi) .* ((MaxFEs-gg)/MaxFEs) + 1 ) .* ones(N,dim);
                flag1 = true;
            else
                sf = 0.5 *( sin(2*pi .* freq(:, ones(1, dim)) .* gg) .* (gg/MaxFEs) + 1 ) .* ones(N,dim);
                flag2 = true;
            end
        else
            ns1 = size(goodF1,1);
            ns1_sum = 0;
            nf1_sum = 0;
            for hh = gg-LP : gg-1
                ns1_sum = ns1_sum + goodF1all(1,hh);
                nf1_sum = nf1_sum + badF1all(1,hh);
            end
            sumS1 = (ns1_sum/(ns1_sum + nf1_sum)) + 0.01;
            ns2 = size(goodF2,1);
            ns2_sum = 0;
            nf2_sum = 0;
            for hh = gg-LP : gg-1
                ns2_sum = ns2_sum + goodF2all(1,hh);
                nf2_sum = nf2_sum + badF2all(1,hh);
            end
            sumS2 = (ns2_sum/(ns2_sum + nf2_sum)) + 0.01;
            p1 = sumS1/(sumS1 + sumS2);
            p2 = sumS2/(sumS2 + sumS1);
            if(p1 > p2)
                sf = 0.5.*( sin(2.*pi.*freq_inti.*gg+pi) .* ((MaxFEs-gg)/MaxFEs) + 1 ) .* ones(N,dim);
                flag1 = true;
            else
                sf = 0.5 *( sin(2*pi .* freq(:, ones(1, dim)) .* gg) .* (gg/MaxFEs) + 1 ) .* ones(N,dim);
                flag2 = true;
            end
        end
    end
    r0 = 1 : N;
    XAll = X;
    [r1, r2] = gnR1R2(N, size(XAll, 1), r0);
    pNP = max(round(p_best_rate * N), 2); 
    randindex = ceil(rand(1, N) .* pNP);
    randindex = max(1, randindex); 
    pbest = X(sorted_index(randindex), :); 
    vi = X + sf(:, ones(1, dim)) .* (pbest - X + X(r1, :) - XAll(r2, :));
    vi = boundConstraint(vi, X, lu);
    J_= mod(floor(rand(N, 1)*dim), dim) + 1;
    J = (J_-1)*N + (1:N)';
    crs = rand(N, dim) < cr(:, ones(1, dim));
    if rand<pb
        best = X(sorted_index(1), :);
        Dis = pdist2(X,best,'euclidean'); % euclidean distance
        [~, idx_ordered] = sort(Dis, 'ascend');
        Neighbour_best_pool = X(idx_ordered(1:SEL), :); %%% including best also so start from 1
        Xsel = Neighbour_best_pool;
        xmean = mean(Xsel);
        C =  1/(SEL-1)*(Xsel - xmean(ones(SEL,1), :))'*(Xsel - xmean(ones(SEL,1), :));
        C = triu(C) + transpose(triu(C,1)); % enforce symmetry
        [R,D] = eig(C);
        if max(diag(D)) > 1e20*min(diag(D))
            tmp = max(diag(D))/1e20 - min(diag(D));
            C = C + tmp*eye(dim);
            [R, D] = eig(C);
        end
        TM = R;
        TM_=R';
        Xr = X*TM;
        vi = vi*TM;
        Ur = Xr;
        Ur(J) = vi(J);
        Ur(crs) = vi(crs);
        ui = Ur*TM_;
    else
        ui = X;
        ui(J) = vi(J);
        ui(crs) = vi(crs);
    end
    children_fitness = [];
    for i=1:size(ui,1)
        children_fitness(i)=fobj(ui(i,:));
        FEs=FEs+1;
    end
    if size(children_fitness,1)==1
       children_fitness = children_fitness';
    end
    for i = 1 : N
        if children_fitness(i) < bestFitness
            bestFitness = children_fitness(i);
            best_pos = ui(i, :);
            best_index = i;
        end
        if FEs > MaxFEs; break; end
    end
    dif = abs(fitness - children_fitness);
    I = (fitness > children_fitness);
    goodCR = cr(I == 1);
    goodF = sf(I == 1);
    goodFreq = freq(I == 1);
    dif_val = dif(I == 1);
    badF = sf(I == 0);
    if flag1 == true
        goodF1 = goodF;
        goodF1all = [goodF1all size(goodF1,1)];
        badF1 = badF;
        badF1all = [badF1all size(badF1,1)];
        goodF2all = [goodF2all 1];
        badF2all = [badF2all 1];
    end
    if flag2 == true
        goodF2 = goodF;
        goodF2all = [goodF2all size(goodF2,1)];
        badF2 = badF;
        badF2all = [badF2all size(badF2,1)];
        goodF1all = [goodF1all 1];
        badF1all = [badF1all 1];
    end
    [fitness, I] = min([fitness, children_fitness], [], 2);
    run_funcvals = [run_funcvals; fitness];
    popold = X;
    popold(I == 2, :) = ui(I == 2, :);
    num_success_params = numel(goodCR);
    if num_success_params > 0
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;
        memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
        if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
            memory_cr(memory_pos)  = -1;
        else
            memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
        end
        if max(goodFreq) == 0 || memory_freq(memory_pos)  == -1
            memory_freq(memory_pos)  = -1;
        else
            memory_freq(memory_pos) = (dif_val' * (goodFreq .^ 2)) / (dif_val' * goodFreq);
        end
        memory_pos = memory_pos + 1;
        if memory_pos > memory_size;  memory_pos = 1; end
    end
    bestFitness = min(fitness);
    Convergence_curve(l)=bestFitness;
    l=l+1;
end 
toc
end

function vi = boundConstraint (vi, pop, lu)

[NP, D] = size(pop);  % the population size and the problem's dimension

%% check the lower bound
xl = repmat(lu(1, :), NP, 1);

pos = vi < xl;
vi(pos) = (pop(pos) + xl(pos)) / 2;

%% check the upper bound
xu = repmat(lu(2, :), NP, 1);
pos = vi > xu;
vi(pos) = (pop(pos) + xu(pos)) / 2;
end

function [r1, r2] = gnR1R2(NP1, NP2, r0)

NP0 = length(r0);

r1 = floor(rand(1, NP0) * NP1) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = (r1 == r0);
    if sum(pos) == 0
        break;
    else % regenerate r1 if it is equal to r0
        r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate r1 in 1000 iterations');
    end
end

r2 = floor(rand(1, NP0) * NP2) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = ((r2 == r1) | (r2 == r0));
    if sum(pos)==0
        break;
    else % regenerate r2 if it is equal to r0 or r1
        r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate r2 in 1000 iterations');
    end
end
end