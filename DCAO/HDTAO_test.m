%% Hierarchical  Dual-Thinking AO algorithm.
function [best_pos, Convergence_curve]=HDTAO(N,MaxFEs,Lb,Ub,dim, fobj)
% Initialization parameters
FEs=0;
it=1;

%% Initialization of the solution set
X=initialization(N,dim,Ub,Lb);
%Calculate the fitness value of the initial solution set
AllFitness = zeros(N,1);
for i=1:N
    AllFitness(i, 1)=fobj(X(i,:));
    FEs=FEs+1;
end
[fmin,~]=min(AllFitness);

%Container
newX=zeros(N,dim);
Convergence_curve=[];
%Record the current optimal solution
bestFitness=fmin;

Lb = Lb * ones(1, dim);
Ub = Ub * ones(1, dim);
newX = zeros(N,dim);
newFitness = zeros(N,1);
eta_qKR = zeros(1,N); %%%%%
golden_ratio = 2/(1 + sqrt(5)); %%%%%
% ngS = max(6,round(N * (golden_ratio/3))); %%%%%
phi_qKR = 0.25 + 0.55 * ((0 + ((1:N)/N)) .^ 0.5); %%%%%


%% Main loop
while FEs<=MaxFEs
    
    ngS = max(6,round(6 + ((FEs/MaxFEs) * 9))); %%%%%
    
    K= 1-((FEs)^(1/6)/(MaxFEs)^(1/6));
    E =1*exp(-4*(FEs/MaxFEs));
    %
    
    % Sort population by fitness values
    [X, AllFitness, ~] = PopSort(X,AllFitness);
    % Reset
    bestInd = 1;
    lamda_t = 0.1 + (0.518 * ((1-(FEs/MaxFEs)^0.5))); %%%%%
    
    for i=1: N
        for j=1:dim
            if rand<K
                newX(i,j) = X(i,j)+E.*X(i,j)*(-1)^FEs;
            else
                newX(i,j)=X(i,j);
            end
        end

        eta_qKR(i) = (round(rand * phi_qKR(i)) + (rand <= phi_qKR(i)))/2; %%%%%
        if i <= ngS %%%%%
            while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= bestInd, break, end, end
            for d = 1:dim
                if rand <= eta_qKR(i)
                    newX(i,d) = X(r1,d) + LnF3(golden_ratio,0.05,1,1); %(Divergent thinking)
                end
            end
        else
            while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= bestInd, break, end, end
            while true, r2 = ngS + round((N - ngS) * rand + 0.5); if r2 ~= i && r2 ~= bestInd && r2 ~= r1, break, end, end
            omega_it = rand; %%%%%
            for d = 1:dim
                if rand <= eta_qKR(i)
                    newX(i,d) = X(bestInd,d) + ((X(r2,d) - X(i,d)) * lamda_t) + ((X(r1,d) - X(i,d)) * omega_it); %(Convergent thinking)
                end
            end
        end
        % Boundary
        newX(i,:) = boundConstraint(newX(i,:),X(i,:),[Lb; Ub]);
        newFitness(i,1) = fobj(newX(i, :));
        FEs = FEs + 1;
        if newFitness(i,1) <= AllFitness(i,1)
            X(i,:) = newX(i,:);
            AllFitness(i,1) = newFitness(i,1);
            if newFitness(i,1) < bestFitness
                bestFitness = newFitness(i,1);
                bestInd = i;
            end
        end
    end
    
    best_pos = X(bestInd,:);
    best_cost = bestFitness;
    Convergence_curve(it) = bestFitness;
    it = it + 1;
end
end

function [sorted_population, sorted_fitness, sorted_index] = PopSort(input_pop,input_fitness)
    [sorted_fitness, sorted_index] = sort(input_fitness,1,'ascend');
    sorted_population = input_pop(sorted_index,:);
end


function Y = LnF3(alpha, sigma, m, n)
    Z = laplacernd(m, n);
    Z = sign(rand(m,n)-0.5) .* Z;
    U = rand(m, n);
    R = sin(0.5*pi*alpha) .* tan(0.5*pi*(1-alpha*U)) - cos(0.5*pi*alpha);
    Y = sigma * Z .* (R) .^ (1/alpha);
end


function x = laplacernd(m, n)
    u1 = rand(m, n);
    u2 = rand(m, n);
    x = log(u1./u2);
end


function vi = boundConstraint(vi, pop, lu)
    % if the boundary constraint is violated, set the value to be the middle
    % of the previous value and the bound
    [NP, ~] = size(pop);  % the population size and the problem's dimension
    % check the lower bound
    xl = repmat(lu(1, :), NP, 1);
    pos = vi < xl;
    vi(pos) = (pop(pos) + xl(pos)) / 2;
    % check the upper bound
    xu = repmat(lu(2, :), NP, 1);
    pos = vi > xu;
    vi(pos) = (pop(pos) + xu(pos)) / 2;
end