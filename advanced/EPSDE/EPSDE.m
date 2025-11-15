%**************************************************************************************************
%Reference:  R. Mallipeddi, P. N. Suganthan, Q. K. Pan, and M. F. Tasgetiren,
%                     " Differential evolution algorithm with ensemble of parameters and
%                     mutation strategies," Applied Soft Computing, in press, 2010.
%
% Note: We obtained the MATLAB source code from the authors, and did some
%           minor revisions in order to solve the 25 benchmark test functions,
%           however, the main body was not changed.
%**************************************************************************************************

function [GBEST, cg_curve]=EPSDE(N,MaxFEs,lb,ub,dim,fobj)
tic;

format long;
format compact;

GBEST = zeros(1,dim);
fes = 0;
cg_curve=[];
t = 1;
lu = [lb .* ones(1, dim); ub .* ones(1, dim)];

% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.

% Define the dimension of the problem
n = dim;

% Main body which was developed by the authors
I_NP = N;
I_D = dim;
FVr_minbound = lu(1, :);
FVr_maxbound = lu(2, :);
Lbound = lu(1, :);
Ubound = lu(2, :);
MAX_FEs = n * 10000;


%-----Initialize population and some arrays-------------------------------

FM_pop = zeros(I_NP, I_D); %initialize FM_pop to gain speed


FM_pop = repmat(FVr_minbound, I_NP, 1) + rand(I_NP, I_D) .* (repmat(FVr_maxbound, I_NP, 1) - repmat(FVr_minbound, I_NP, 1));


FM_popold  = zeros(size(FM_pop)); % toggle population
val   = zeros(I_NP, 1);    %create and reset the "cost array"
FVr_bestmem = zeros(1, I_D); % best population member ever
I_nfeval  = 0;      % number of function evaluations

%------Evaluate the best member after initialization----------------------

I_best_index = 1;     % start with first population member
val(1)   = fobj(FM_pop(I_best_index, :));
fes = fes + 1;
F_bestval = val(1);     % best objective function value so far
I_nfeval = I_nfeval + 1;
for k = 2 : I_NP       % check the remaining members
    val(k) = fobj(FM_pop(k, :));
    fes = fes + 1;
    I_nfeval = I_nfeval + 1;
    if (val(k) < F_bestval)

        I_best_index = k;    % save its location
        F_bestval  = val(k);
    end
end
FVr_bestmem = FM_pop(I_best_index, :); % best member of current iteration

%------DE-Minimization---------------------------------------------
%------FM_popold is the population which has to compete. It is--------
%------static through one iteration. FM_pop is the newly--------------
%------emerging population.----------------------------------------
I_iter = 1;

FF = [0.4; 0.5; 0.6; 0.7; 0.8; 0.9];
CR = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9];

spara = randperm(3)';
Para = zeros(I_NP, 3);
inde = 1 : I_NP;
inde = inde';
RR = zeros(I_NP, 1);
PPara = [];

% while (I_nfeval < MAX_FEs)
while ( fes <= MaxFEs )

%     'EPSDE'
%     t

    if(I_iter == 1)

%         Para(inde, :) = [spara(randi(size(inde, 1), 1, [1, size(spara, 1)])), CR(randi(size(inde, 1), 1, [1, size(CR, 1)])), FF(randi(size(inde, 1), 1, [1, size(FF, 1)]))];
        Para(inde, :) = [spara(randi(size(spara, 1), 1, size(inde, 1))), CR(randi(size(CR, 1), 1, size(inde, 1))), FF(randi(size(FF, 1), 1, size(inde, 1)))];

             
    else
        for k = 1 : size(inde, 1)

            if(rand <= RATE && ~isempty(PPara))

%                 RR(k) = randi(1, 1, [1, size(PPara, 1)]);
                RR(k) = randi(size(PPara, 1), 1, 1);
%                 Para(inde(k), :) = PPara(randi(1, 1, [1, size(PPara, 1)]), :);
                Para(inde(k), :) = PPara(randi(size(PPara, 1)), :);


            else

                RR(k) = 0;
%                 Para(inde(k), :) = [spara(randi(1, 1, [1, size(spara, 1)])), CR(randi(1, 1, [1, size(CR, 1)])), FF(randi(1, 1, [1, size(FF, 1)]))];
                Para(inde(k), :) = [spara(randi(size(spara, 1), 1)), CR(randi(size(CR, 1), 1)), FF(randi(size(FF, 1), 1))];

            end
        end
    end

    RRR = [];
    count = 0;
    FM_popold = FM_pop;

    for i = 1 : I_NP

        FM_mui = rand(1, I_D) < Para(i, 2);
        dd = find(FM_mui == 1);
        if isempty(dd)
            ddd = ceil(rand * I_D);
            FM_mui(ddd) = 1;
        end
        FM_mpo = FM_mui < 0.5;
        FM_bm = FVr_bestmem;
        para(i, :) = normrnd(Para(i, 3), 0.001, 1, I_D);
        if(Para(i, 1) == 1)
            %DE/best/2/bin
            ind = randperm(I_NP);
            FM_pm3 = FM_popold(ind(1), :);
            FM_pm4 = FM_popold(ind(2), :);
            FM_pm5 = FM_popold(ind(3), :);
            FM_pm6 = FM_popold(ind(4), :);
            FM_ui(i, :) = FM_bm +(FM_pm3 - FM_pm4 + FM_pm5 - FM_pm6) .* para(i, :);
            FM_ui(i, :) = FM_popold(i, :) .* FM_mpo + FM_ui(i, :) .* FM_mui;
        end
        if(Para(i, 1) == 2)
            %DE/rand/1/bin
            ind = randperm(I_NP);
            FM_pm7 = FM_popold(ind(1), :);
            FM_pm8 = FM_popold(ind(2), :);
            FM_pm9 = FM_popold(ind(3), :);
            FM_ui(i, :) = FM_pm7 +para(i, :) .* (FM_pm8 - FM_pm9);
            FM_ui(i, :) = FM_popold(i, :) .* FM_mpo+ FM_ui(i, :) .* FM_mui;  % crossover
        end
        if(Para(i, 1) == 3)
            % DE/current-to-rand/1/bin/
            ind = randperm(I_NP);
            FM_pm21 = FM_popold(ind(1), :);
            FM_pm22 = FM_popold(ind(2), :);
            FM_pm23 = FM_popold(ind(3), :);
            FM_ui(i, :) = FM_popold(i, :) + rand(1, I_D) .* (FM_pm21-FM_popold(i, :)) + para(i, :) .* (FM_pm22 - FM_pm23);  % differential variation
        end

        if(FM_ui(i, :) < Lbound)
            FM_ui(i, :) = FVr_minbound + (FVr_maxbound-FVr_minbound) .* rand(1, I_D);

        end
        if(FM_ui(i, :) > Ubound)
            FM_ui(i, :) = FVr_minbound + (FVr_maxbound-FVr_minbound) .* rand(1, I_D);
        end

    end

    tempval = zeros(N, 1);
    for i  = 1 : N
        tempval(i) = fobj(FM_ui(i, :));
        fes = fes + 1;
    end

    for i = 1 : I_NP

        I_nfeval = I_nfeval + 1;
        if (tempval(i) < val(i))
            FM_pop(i, :) = FM_ui(i, :);
            val(i) = tempval(i);
            PPara = [Para(i, :); PPara];
            if(RR(i) ~= 0)
                RRR = [RRR; RR(i)];
            end
            if (tempval(i) < F_bestval)
                F_bestval = tempval(i);
                FVr_bestmem = FM_ui(i, :);
                I_best_index = i;
            end

        else

            count = count + 1;

        end

    end

    PPara(RRR, :) = [];
    rate(I_iter, 1) = count / I_NP;
    if(I_iter > 10)
        RATE = mean(rate((I_iter-10) : I_iter), 1);
    else
        RATE = mean(rate, 1);
    end

    I_iter = I_iter + 1;

    cg_curve(t) = F_bestval;
    t = t + 1;
end

toc;
