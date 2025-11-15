function [Best_pos,Convergence_curve]=DWPLO(N,MaxFEs,lb,ub,dim,fobj)

%% Initialization
FEs = 0;
it = 1;
fitness=inf*ones(N,1);
fitness_new=inf*ones(N,1);

X=initialization(N,dim,ub,lb);
V=ones(N,dim);
X_new=zeros(N,dim);

for i=1:N
    fitness(i)=fobj(X(i,:));
    FEs=FEs+1;
end

[fitness, SortOrder]=sort(fitness);
X=X(SortOrder,:);
Bestpos=X(1,:);
Bestscore=fitness(1);

Convergence_curve=[];
Convergence_curve(it)=Bestscore;

%% Main loop
while FEs <= MaxFEs
    
    X_sum=sum(X,1);
    X_mean=X_sum/N;
    w1=tansig((FEs/MaxFEs)^4);
    w2=exp(-(2*FEs/MaxFEs)^3);

    F = 0.5;  % Differential evolution scaling factor
    CR = 0.5;  % Crossover probability
    
    for i=1:N
        
        a=rand()/2+1;
        V(i,:)=1*exp((1-a)/100*FEs);
        LS=V(i,:);

        GS=Levy(dim).*(X_mean-X(i,:)+(lb+rand(1,dim)*(ub-lb))/2);
        if rand < CR
            r1 = randi([1 N]);  % Randomly select r1
            r2 = randi([1 N]);  % Randomly select r2
            r3 = randi([1 N]);  % Randomly select r3
            while r1 == i, r1 = randi([1 N]); end
            while r2 == i || r2 == r1, r2 = randi([1 N]); end
            while r3 == i || r3 == r1 || r2 == r3, r3 = randi([1 N]); end
            X_new(i, :) = X(r1, :) + F * (X(r2, :) - X(r3, :));
        else
            X_new(i, :) = X(i, :) + (w1 * LS + w2 * GS) .* rand(1, dim);
        end
    end
    
    E =sqrt(FEs/MaxFEs);
    A=randperm(N);
    for i=1:N
        for j=1:dim
            if (rand<0.05) && (rand<E)
                X_new(i,j)=X(i,j)+sin(rand*pi)*(X(i,j)-X(A(i),j));
            elseif rand > (0.5*E + 0.5)
                X_new(i,j)= WeilbullFlight(X(i,j),Bestpos(j));
            end
        end
        Flag4ub=X_new(i,:)>ub;
        Flag4lb=X_new(i,:)<lb;
        X_new(i,:)=(X_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        fitness_new(i)=fobj(X_new(i,:));
        FEs=FEs+1;
        if fitness_new(i)<fitness(i)
            X(i,:)=X_new(i,:);
            fitness(i)=fitness_new(i);
        end
    end
    [fitness, SortOrder]=sort(fitness);
    X=X(SortOrder,:);
    if fitness(1)<Bestscore
        Bestpos=X(1,:);
        Bestscore=fitness(1);
    end
    it = it + 1;
    Convergence_curve(it)=Bestscore;
    Best_pos=Bestpos;
end

end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);
step=u./abs(v).^(1/beta);
o=step;
end

function Xnew=WeilbullFlight(X,Best_agent)
if rand <=0.25
    X=X+sign(rand()-0.5).*wblrnd(1,1);
else
    if  isequal (Best_agent,X)
        step=0.1*sign(rand()-0.5);
        X=X+step.*wblrnd(0.5,1);
    else
        step=0.5*sign(rand()-0.5)*norm(Best_agent-X);
        X=X+step.*wblrnd(0.5,1); %0.5
    end
end
Xnew = X;
end