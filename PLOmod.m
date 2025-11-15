function [best_pos,Convergence_curve]=PLOmod(N,MaxFEs,lb,ub,dim,fobj)

%% Initialization
FEs = 0;
it = 1;
AllFitness=inf*ones(N,1);
newAllFitness=inf*ones(N,1);

X=initialization(N,dim,ub,lb);
V=ones(N,dim);
newX=zeros(N,dim);

for i=1:N
    AllFitness(i)=fobj(X(i,:));
    FEs=FEs+1;
end

[AllFitness, SortOrder]=sort(AllFitness);
X=X(SortOrder,:);
Best_pos=X(1,:);
bestFitness=AllFitness(1);

rec = 1;
rec_num = 10;
rM = zeros(N,dim,rec_num); %record history positions
rM_cos = zeros(1,N,rec_num);

Convergence_curve=[];
Convergence_curve(it)=bestFitness;

%% Main loop
while FEs <= MaxFEs
    
    X_sum=sum(X,1);
    X_mean=X_sum/N;
    w1=tansig((FEs/MaxFEs)^4);
    w2=exp(-(2*FEs/MaxFEs)^3);
    
    if rec == 1 %record the first generation of positions
        rM(:,:,rec) = X;
        rM_cos(1,:,rec) = AllFitness;
        rec = rec + 1;
    end

    for i=1:N
        
        a=rand()/2+1;
        V(i,:)=1*exp((1-a)/100*FEs);
        LS=V(i,:);

        GS=Levy(dim).*(X_mean-X(i,:)+(lb+rand(1,dim)*(ub-lb))/2);
        newX(i, :) = X(i, :) + (w1 * LS + w2 * GS) .* rand(1, dim); %%1
    end
    
    E =sqrt(FEs/MaxFEs);
    A=randperm(N);
    for i=1:N
        for j=1:dim
            if (rand<0.05) && (rand<E)
                newX(i,j)=X(i,j)+sin(rand*pi)*(X(i,j)-X(A(i),j)); %%2
            end
        end

%% -------------------Migration behavior-------------------%
        r = rand();
        m=2*sin(r+pi/2);
        s = randi([1,30],1);
        tempFitness=AllFitness(s);
        ori_value = rand(1,dim);cauchy_value = tan((ori_value-0.5)*pi);
        if AllFitness(i) < tempFitness
            newXX(i,:)=X(i,:)+cauchy_value(:,dim).* (X(i,:)-Best_pos);%%3
        else
            newXX(i,:)=X(i,:)+cauchy_value(:,dim).* (Best_pos-m.*X(i,:));%%4
            % newXX(i, :) = X(i, :) + (w1 * LS + w2 * GS) .* rand(1, dim); %%1
        end


        Flag4ub=newX(i,:)>ub;
        Flag4lb=newX(i,:)<lb;
        newX(i,:)=(newX(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        newAllFitness(i)=fobj(newX(i,:));
        FEs=FEs+1;

        Flag4ub=newXX(i,:)>ub;
        Flag4lb=newXX(i,:)<lb;
        newX(i,:)=(newXX(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        tempFitness=fobj(newXX(i,:));
        FEs=FEs+1;

        if tempFitness<newAllFitness(i)
            newX(i,:)=newXX(i,:);
            newAllFitness(i)=tempFitness;
        end

        % Cryptobiosis mechanism
        rM(i,:,rec) = newX(i,:);
        rM_cos(1,i,rec) = newAllFitness(i);

        if newAllFitness(i)<AllFitness(i)
            X(i,:)=newX(i,:);
            AllFitness(i)=newAllFitness(i);
        end
    end
    [AllFitness, SortOrder]=sort(AllFitness);
    X=X(SortOrder,:);
    if AllFitness(1)<bestFitness
        Best_pos=X(1,:);
        bestFitness=AllFitness(1);
    end

    rec = rec + 1;
    % Cryptobiosis mechanism
    if rec > rec_num || FEs>=MaxFEs
        [lcost,Iindex] = min(rM_cos, [] ,3);
        for i=1:N
            X(i,:) = rM(i,:,Iindex(i));
        end
        AllFitness = lcost;
        rec = 1;
    end

    it = it + 1;
    Convergence_curve(it)=bestFitness;
    best_pos=Best_pos;
end

end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);
step=u./abs(v).^(1/beta);
o=step;
end
