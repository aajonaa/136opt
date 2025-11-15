
function [best_Location,CNVG]=yingwu(N,Max_FEs,lb,ub,dim,fobj)

% disp('鹦鹉 is now tackling your problem，芝士鹦鹉')
tic
FEs = 0;
% initialize the location and Energy of the rabbit
best_Location=zeros(1,dim);
best_fitness=inf;

%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);
X_fitness=zeros(1,N);
X_new=X;
fitness_new=zeros(1,N);

for i=1:N
    FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
    X_fitness(i)=fobj(X(i,:));
    FEs = FEs + 1;

    if X_fitness(i)<best_fitness
        best_fitness=X_fitness(i);
        best_Location=X(i,:);
    end
end

CNVG=[];

iter = 1;

while FEs < Max_FEs
    alpha=rand()/5;
    sita=rand()/pi;
    for j = 1:N
        St=randi([1,4]);

        if St==1
            X_new(j, :) = (X(j, :) - best_Location) .* Levy(dim) + rand(1) * mean(X(j, :)) * (1 - i / Max_FEs) ^ (2 * i / Max_FEs);
        elseif St==2
            X_new(j, :) = X(j, :) + best_Location .* Levy(dim) + randn() * (1 - i / Max_FEs) * ones(1, dim);
        elseif St==3
            H = rand(1);
                if H < 0.5
                    X_new(j, :) = X(j, :) + alpha * (1 - i / Max_FEs) * (X(j, :) - mean(X(j, :)));
                else
                    X_new(j, :) = X(j, :) + alpha * (1 - i / Max_FEs) * exp(-j / (rand(1) * Max_FEs));
                end
        else
             X_new(j, :) = X(j, :) + rand() * cos((pi *i )/ (2 * Max_FEs)) * (best_Location - X(j, :)) - cos(sita) * (i / Max_FEs) ^ (2 / Max_FEs) * (X(j, :) - best_Location);
        end
        FU=X_new(i,:)>ub;
        FL=X_new(i,:)<lb;
        X_new(i,:)=(X_new(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        fitness_new(i)=fobj(X_new(i,:));
        FEs = FEs + 1;
        
        if fitness_new(i)<best_fitness
            best_fitness=fitness_new(i);
            best_Location=X_new(i,:);
        end
    end
    X = X_new;
    fitness = fitness_new;

    CNVG(iter)=best_fitness;
    iter = iter + 1;

end
toc
end

function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) *sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1 / beta);
    o = step;
end