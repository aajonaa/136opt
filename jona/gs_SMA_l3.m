% Slime Mold Algorithm modified by Jona 2023-11-16.
% 1. gamma = 2. 2. add rand() * 1 / N. 3.shift the value 4. golden sine
% mechanisim
function [best_pos,convergence_curve]=gs_SMA_l3(N,Max_FEs,lb,ub,dim,fobj)
tic
disp('gs_SMA_l3 is now tackling your problem')

% initialize position
best_pos=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(1,N);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
% X=initialization(N,dim,ub,lb);

% %% sine map
% for i = 1:N
%     if i == 1
%         for j = 1:dim
%             X(i, j) = (ub - lb) * rand + lb;
%         end
%     else
%         a = rand();
%         for j = 1:dim
%             X(i, j) = (4 * sin(pi * X(i - 1, j))) / (4 * a);
%         end
%     end
% end

% %% logistic tent map
% % mu = 4;
% % gamma = 0.5;
% for i = 1:N
%     mu = 4 * rand();
%     gamma = 0.5 * rand();
%     for j = 1:dim
%         x = lb + (ub - lb) * rand();
% %         X(i, j) = mu * x * (1 - x) * gamma * min(x, 1-x);
%         X(i, j) = mu * (x - lb) / (ub - lb) * (1 - (x - lb) / (ub - lb)) * gamma * min((x - lb) / (ub - lb), 1 - (x - lb) / (ub - lb));
%     end
% end

%% logistic map
% mu = 4;
for i = 1:N
    if i == 1
        for j = 1:dim
            X(i, j) = (ub - lb) * rand() + lb;
        end
    else
        mu = 4 * rand();
        for j = 1:dim
            X(i, j) = mu * (X(i-1, j) - lb) / (ub - lb) * (1 - (X(i-1, j) - lb) / (ub - lb)) * (ub - lb) + lb;
        end
    end
end

% %% tent map
% % gamma = 0.5;
% for i = 1:N
%     if i == 1
%         for j = 1:dim
%             X(i, j) = (ub - lb) * rand() +lb;
%         end
%     else
% %         gamma = 0.5 * rand();
%         for j = 1:dim
%             X(i, j) = 2 * min((X(i - 1, j) - lb) / (ub - lb), 1 - (X(i-1, j) - lb) / (ub - lb)) * (ub - lb) + lb + rand() * 1 / N;
%         end
%     end
% end

% %% logistic tent map
% % mu = 4;
% % gamma = 0.5;
% for i = 1:N
%     if i == 1
%         for j = 1:dim
%             X(i, j) = (ub - lb) * rand() + lb;
%         end
%     else
%         mu = 4 * rand();
%         gamma = 0.5 * rand();
%         for j = 1:dim
% %         X(i, j) = mu * x * (1 - x) * gamma * min(x, 1-x);
%             X(i, j) = mu * (X(i-1, j) - lb) / (ub - lb) * (1 - (X(i - 1, j) - lb) / (ub - lb)) * gamma * min((X(i-1, j) - lb) / (ub - lb), 1 - (X(i-1, j) - lb) / (ub - lb));
%         end
%     end
% end
% %%

convergence_curve=[];
it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

% Main loop
FEs = 0;



a = pi;
b = -pi;
tao = (sqrt(5) - 1) / 2;
c1 = a * tao + b * (1 - tao);
c2 = a * (1 - tao) + b * tao;




while  FEs <= Max_FEs
    
    %sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
        FEs = FEs + 1;
    end
    
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);

    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero

    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        best_pos=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end


    %% golden sine
%     a = pi;
%     b = -pi;
%     tao = (sqrt(5) - 1) / 2;
%     c1 = a * tao + b * (1 - tao);
%     c2 = a * (1 - tao) + b * tao;
%     x_new = X(q, :) .* abs(sin(r1)) - r2 * sin(r1) * abs(c1 .* best_pos - c2 .* X(q, :));
%     for i = 1:N
%         r1 = 2 * pi * rand();
%         r2 =  pi * rand();
%         for j = 1:dim
%             X(i, j) = X(i, j) * abs(sin(r1)) - r2 * sin(r1) * abs(c1 * best_pos(j) - c2 * X(i, j));
% %             X1(i, j) = best_pos(j) * abs(sin(r1)) - r2 * sin(r1) * abs(c1 * best_pos(j) - c2 * best_pos(j));
%         end
% %         FEs = FEs + 1;
% %         fitness(i) = fobj(X1(i, :));
%     end
% 
% 
%     for i=1:size(X,1)
%         
%         % Check if solutions go outside the search spaceand bring them back
%         Flag4ub=X(i,:)>ub;
%         Flag4lb=X(i,:)<lb;
%         X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
%         
%         % Calculate the objective values
%         FEs = FEs + 1;
%         Destination_values(1,i)=fobj(X(i,:));
%         
%         % Update the destination if there is a better solution
%         if Destination_values(1,i)<Destination_fitness
%             best_pos=X(i,:);
%             Destination_fitness=Destination_values(1,i);
%             b = c2;
%             c2 = c1;
%             c1 = a * tao + b * (1 - tao);
%         else
%             a = c1;
%             c1 = c2;
%             c2 = a * (1 - tao) + b * tao;
%         end
%         if c1 == c2
%             a = pi * rand();
%             b = -pi * rand();
%             c1 = a * tao + b * (1 - tao);
%             c2 = a * (1 - tao) + b * tao;
%         end
%         
%     end

%     ind = fitness < AllFitness;
%     disp(ind);
%     X(ind, :) = X1(ind, :);
%     [sorted_fitness, index] = sort(fitness);
%     if sorted_fitness(1) < Destination_fitness
%         best_pos = X1(index(1), :);
%         b = c2;
%         c2 = c1;
%         c1 = a * tao + b * (1 - tao);
%     else
%         a = c1;
%         c1 = c2;
%         c2 = a * (1 - tao) + b * tao;
%     end
%     if c1 == c2
%         a = pi * rand();
%         b = -pi * rand();
%         c1 = a * tao + b * (1 - tao);
%         c2 = a * (1 - tao) + b * tao;
%     end
    


%         if fitness < Destination_fitness
%             Destination_fitness = fitness;
%             best_pos = X1(i, :);
%             b = c2;
%             c2 = c1;
%             c1 = a * tao + b * (1 - tao);
%         else
%             a = c1;
%             c1 = c2;
%             c2 = a * (1 - tao) + b * tao;
%         end
%         if c1 == c2
%             a = pi * rand();
%             b = -pi * rand();
%             c1 = a * tao + b * (1 - tao);
%             c2 = a * (1 - tao) + b * tao;
%         end
%     end


%     FEs = FEs + 1;
%     if fobj(x_new) < bestFitness
%         best_pos = x_new;
%         b = c2;
%         c2 = c1;
%         c1 = a * tao + b * (1 - tao);
%     else
%         b = c1;
%         c1 = c2;
%         c2 = a * (1 - tao) + b * tao;
%     end


    %%

%     if rand < ((1-FEs) / Max_FEs)^(1 - tan(pi * (rand - 0.5)) * s / Max_FEs)
%         m = 1;
%     else
%         m = 0;
%     end

    
    a = atanh(-(FEs/Max_FEs)+1);   %Eq.(2.4)
    b = 1-FEs/Max_FEs;
    % Update the Position of search agents
    for i=1:N
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(2.2)
            vb = unifrnd(-a,a,1,dim);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,dim);

            r1 = 2 * pi * rand();
            r2 =  pi * rand();
            s = 1;

            for j=1:dim
                r = rand();
                A = randi([1,N]);  % two positions randomly selected from population
                B = randi([1,N]);
                if r<p    %Eq.(2.1)
                    X(i,j) = best_pos(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    if rand < ((1-FEs) / Max_FEs)^(1 - tan(pi * (rand - 0.5)) * s / Max_FEs)
                        m = 1;
                        X(i,j) = vc(j)*X(i,j);
                        a = c1;
                        c1 = c2;
                        c2 = a * (1 - tao) + b * tao;
                    else
                        m = 0;
                        X(i, j) = X(i, j) * abs(sin(r1)) - r2 * sin(r1) * abs(c1 * best_pos(j) - c2 * X(i, j));
                        b = c2;
                        c2 = c1;
                        c1 = a * tao + b * (1 - tao);
                    end
%                     X(i,j) = vc(j)*X(i,j);
%                     X(i, j) = X(i, j) * abs(sin(r1)) - r2 * sin(r1) * abs(c1 * best_pos(j) - c2 * X(i, j));
%                     b = c2;
%                     c2 = c1;
%                     c1 = a * tao + b * (1 - tao);
                    if c1 == c2
                        a = pi * rand();
                        b = -pi * rand();
                        c1 = a * tao + b * (1 - tao);
                        c2 = a * (1 - tao) + b * tao;
                    end
                end
            end
        end
    end
    convergence_curve(it)=Destination_fitness;
    it=it+1;
end
toc
end
