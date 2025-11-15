function [best_pos,Convergence_curve]=PLO(N,MaxFEs,lb,ub,dim,fobj)

%% Initialization
FEs = 0;
it = 1;
AllFitness=inf*ones(N,1);
newFitness=inf*ones(N,1);

X=initialization(N,dim,ub,lb);
V=ones(N,dim);
newX=zeros(N,dim);

for i=1:N
    AllFitness(i)=fobj(X(i,:));
    FEs=FEs+1;
end

[AllFitness, SortOrder]=sort(AllFitness);
X=X(SortOrder,:);
Bestpos=X(1,:);
bestFitness=AllFitness(1);

Convergence_curve=[];
Convergence_curve(it)=bestFitness;

%% Main loop
while FEs <= MaxFEs
    
    X_sum=sum(X,1);
    X_mean=X_sum/N;
    w1=tansig((FEs/MaxFEs)^4);
    w2=exp(-(2*FEs/MaxFEs)^3);

    % Aurora oval walk
    for i=1:N
        
        a=rand()/2+1;
        V(i,:)=1*exp((1-a)/100*FEs);
        LS=V(i,:);

        GS=Levy(dim).*(X_mean-X(i,:)+(lb+rand(1,dim)*(ub-lb))/2);
        newX(i, :) = X(i, :) + (w1 * LS + w2 * GS) .* rand(1, dim); %%1
    end
    
    % Particle collision
    E =sqrt(FEs/MaxFEs);
    A=randperm(N);
    for i=1:N
        for j=1:dim
            if (rand<0.05) && (rand<E)
                newX(i,j)=X(i,j)+sin(rand*pi)*(X(i,j)-X(A(i),j)); %%2
            end
        end
        Flag4ub=newX(i,:)>ub;
        Flag4lb=newX(i,:)<lb;
        newX(i,:)=(newX(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        newFitness(i)=fobj(newX(i,:));
        FEs=FEs+1;

        if newFitness(i)<AllFitness(i)
            X(i,:)=newX(i,:);
            AllFitness(i)=newFitness(i);
        end
    end
    [AllFitness, SortOrder]=sort(AllFitness);
    X=X(SortOrder,:);
    if AllFitness(1)<bestFitness
        Bestpos=X(1,:);
        bestFitness=AllFitness(1);
    end

    it = it + 1;
    Convergence_curve(it)=bestFitness;
    best_pos=Bestpos;
end

end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);
step=u./abs(v).^(1/beta);
o=step;
end