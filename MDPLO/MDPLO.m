function [best_pos,Convergence_curve]=MDPLO(N,MaxFEs,lb,ub,dim,fobj)
%% Add migration strategy inspired by BKA algorithm and Divergent thinking strategy inspired by DCS algorithm.

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

Convergence_curve=[];
Convergence_curve(it)=bestFitness;

phi_qKR = 0.25 + 0.55 * ((0 + ((1:N)/N)) .^ 0.5); %%%%%
%% Main loop
while FEs <= MaxFEs
    [X, AllFitness, ~] = PopSort(X,AllFitness);
    
    X_sum=sum(X,1);
    X_mean=X_sum/N;
    w1=tansig((FEs/MaxFEs)^4);
    w2=exp(-(2*FEs/MaxFEs)^3);


    for i=1:N

        a=rand()/2+1;
        V(i,:)=1*exp((1-a)/100*FEs);
        LS=V(i,:);

        GS=Levy(dim).*(X_mean-X(i,:)+(lb+rand(1,dim)*(ub-lb))/2);


        %% -------------------Migration behavior-------------------%
        r = rand();
        m=2*sin(r+pi/2);
        s = randi([1,30],1);
        tempFitness=AllFitness(s);
        ori_value = rand(1,dim);cauchy_value = tan((ori_value-0.5)*pi);
        if AllFitness(i) < tempFitness
            newX(i,:) = X(i, :) + (X(i,:)-Best_pos) .* m;%%3 changed the cauchy_value as m
        else
            newX(i, :) = X(i, :) + (w1 * LS + w2 * GS) .* ori_value; %%1
        end

        
    end
    
    A=randperm(N);
    [~, ind] = sort(AllFitness);
    bestInd = ind(1);
    
    %% -------------------Divergent thinking-------------------%
    lamda_t = 0.1 + (0.518 * ((1-(FEs/MaxFEs)^0.5))); %%%%%
    for i=1:N
        jrand = floor(dim * rand + 1); %%%%%
        eta_qKR(i) = (round(rand * phi_qKR(i)) + (rand <= phi_qKR(i)))/2; %%%%%
        omega_it = rand();
        while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= bestInd, break, end, end
        while true, r2 = 6 + round((N - 6) * rand + 0.5); if r2 ~= i && r2 ~= bestInd && r2 ~= r1, break, end, end
        for j=1:dim
            if rand <= eta_qKR(i) || j == jrand
                newX(i,j) = Best_pos(j) + ((X(r2,j) - X(i,j)) * lamda_t) + ((X(r1,j) - X(i,j)) * omega_it);
            else
                newX(i,j)=X(i,j)+sin(rand*pi)*(X(i,j)-X(A(i),j)); %%2
            end
        end
        
        
        %% -------------------The Migration-------------------%
        ori_value = rand(1,dim);cauchy_value = tan((ori_value-0.5)*pi);
        newXX(i,:)=X(i,:)+cauchy_value(:,dim).* (Best_pos-m.*X(i,:));%%4
        
        Flag4ub=newXX(i,:)>ub;
        Flag4lb=newXX(i,:)<lb;
        newXX(i,:)=(newXX(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        FEs=FEs+1;


        Flag4ub=newX(i,:)>ub;
        Flag4lb=newX(i,:)<lb;
        newX(i,:)=(newX(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        newAllFitness(i)=fobj(newX(i,:));
        FEs=FEs+1;


        if fobj(newXX(i, :))<newAllFitness(i)
            newX(i,:)=newXX(i,:);
            newAllFitness(i)=fobj(newXX(i, :));
        end
        
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

function [sorted_population, sorted_fitness, sorted_index] = PopSort(input_pop,input_fitness)
    [sorted_fitness, sorted_index] = sort(input_fitness,1,'ascend');
    sorted_population = input_pop(sorted_index,:);
end
