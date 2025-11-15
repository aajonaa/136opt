function [best_pos, Convergence_curve]=AO5(N,MaxFEs,Lb,Ub,dim, fobj)
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
[fmin,ind]=min(AllFitness);

%Container
newX=zeros(N,dim);
Fitnorm=zeros(1,N);
Convergence_curve=[];
%Record the current optimal solution
best=X(ind,:);
bestFitness=fmin;

pc = 0.5; %%%%%
Lb = Lb * ones(1, dim);
Ub = Ub * ones(1, dim);
nextX = zeros(N,dim);
newX = zeros(N,dim);
newFitness = zeros(N,1);
eta_qKR = zeros(1,N); %%%%%
% Golden ratio
golden_ratio = 2/(1 + sqrt(5)); %%%%%
% High-performing individuals
ngS = max(6,round(N * (golden_ratio/3))); %%%%%
% Ranking-based self-improvement
phi_qKR = 0.25 + 0.55 * ((0 + ((1:N)/N)) .^ 0.5); %%%%%


%% Main loop
while FEs<=MaxFEs
    
    K= 1-((FEs)^(1/6)/(MaxFEs)^(1/6));
    E =1*exp(-4*(FEs/MaxFEs));
    %
    
    % Sort population by fitness values
    [X, AllFitness, ~] = PopSort(X,AllFitness);
    % Reset
    bestInd = 1;
    % Compute social impact factor
    lamda_t = 0.1 + (0.518 * ((1-(FEs/MaxFEs)^0.5))); %%%%%
    
    for i=1: N
        Fitnorm(i)= (AllFitness(i)-min(AllFitness))/(max(AllFitness)-min(AllFitness));
        for j=1:dim
            if rand<K
%                 if rand<0.5
                    newX(i,j) = X(i,j)+E.*X(i,j)*(-1)^FEs;
%                 else
%                     newX(i,j) = X(i,j)+E.*best(j)*(-1)^FEs;
%                 end
            else
                newX(i,j)=X(i,j);
            end
%             if rand<Fitnorm(i)
%                 A=randperm(N);
%                 beta=(rand/2)+0.1;
%                 newX(i,j)=X(A(3),j)+beta.*(X(A(1),j)-X(A(2),j)); 
%             end
        end
%         
%         newX(i,:)=Mutation(newX(i,:),X(i,:),best,dim);
%         newX(i,:)=Transborder_reset(newX(i,:),Ub,Lb,dim,best);
%         
%         tFitness=fobj(newX(i,:));
%         FEs=FEs+1;
%         if tFitness<AllFitness(i)
%             X(i,:)= newX(i,:);
%             AllFitness(i)=tFitness;
%         end


        % Compute differentiated knowledge-acquisition rate
        eta_qKR(i) = (round(rand * phi_qKR(i)) + (rand <= phi_qKR(i)))/2; %%%%%
        jrand = floor(dim * rand + 1); %%%%%
        nextX(i,:) = X(i,:);
%         if i == N && rand < pc
            % Low-performing
%             nextX(i,:) = Lb + rand * (Ub - Lb);
%         elseif i <= ngS %%%%%
%         if i <= ngS %%%%%
            % High-performing
%             while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= bestInd, break, end, end
%             for d = 1:dim
%                 if rand <= eta_qKR(i) || d == jrand
%                     nextX(i,d) = X(r1,d) + LnF3(golden_ratio,0.05,1,1);
%                 end
%             end
%         else
            % Average-performing
            while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= bestInd, break, end, end
            while true, r2 = ngS + round((N - ngS) * rand + 0.5); if r2 ~= i && r2 ~= bestInd && r2 ~= r1, break, end, end
            % Compute learning ability
            omega_it = rand; %%%%%
            for d = 1:dim
                if rand <= eta_qKR(i) || d == jrand
                    nextX(i,d) = X(bestInd,d) + ((X(r2,d) - X(i,d)) * lamda_t) + ((X(r1,d) - X(i,d)) * omega_it);
                end
            end
%         end
        % Boundary
        nextX(i,:) = boundConstraint(nextX(i,:),X(i,:),[Lb; Ub]);
        newX(i, :) = nextX(i, :);
        newFitness(i,1) = fobj(nextX(i, :));
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
    
%     [fmin,ind]=min(AllFitness);
%     if fmin<bestFitness
%         best=X(ind,:);
%         bestFitness=fmin;
%     end
%     
%     Convergence_curve(it)=bestFitness;
%     best_pos = best;
%     it=it+1;
end
end

function z=Mutation(z,x,b,dim)
for j=1:dim
    if rand<0.05
        z(j)=x(j);
    end
    if rand<0.2
        z(j)=b(j);
    end
end
end

function z=Transborder_reset(z,ub,lb,dim,best)
for j=1:dim
    if z(j)>ub || z(j)<lb
        
        z(j)=best(j);
        
    end
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
    %
    % Version: 1.1   Date: 11/20/2007
    % Written by Jingqiao Zhang, jingqiao@gmail.com
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