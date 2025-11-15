%%

function bestFitness=AOPara(N,Max_FEs,lb,ub,dim,fobj, k)

FEs=0;
it=1;
Fitnorm=zeros(1,N);
% Convergence_curve=[];
Fitness = ones(1, N) * inf;

X=initialization(N,dim,ub,lb);

for i=1:N
    Fitness(i)=fobj(X(i,:));
    FEs=FEs+1;
end

[fit_min,minfit_idx]=min(Fitness);
%
newX=zeros(N,dim);
best_position=X(minfit_idx,:);
bestFitness=fit_min;
constant = k;

while FEs<=Max_FEs
    
    K= 1-((FEs)^constant/(Max_FEs)^constant);
    E =1*exp(-4*(FEs/Max_FEs));

    for i=1: N
        
        Fitnorm(i)= (Fitness(i)-min(Fitness))/(max(Fitness)-min(Fitness));
        for j=1:dim
 
            if rand<K
                if rand<0.5
                    newX(i,j) = X(i,j)+E * X(i,j)*(-1)^FEs;
                else
                    newX(i,j) = X(i,j)+E * best_position(j)*(-1)^FEs;
                end
            else
                newX(i,j)=X(i,j);
            end
            
            if rand<Fitnorm(i)
                A=randperm(N);
                beta=(rand/2)+0.1;
                newX(i,j)=X(A(3),j)+beta *(X(A(1),j)-X(A(2),j));
                
            end
        end
        
        newX(i,:)=Mutation(newX(i,:),X(i,:),best_position,dim);
        newX(i,:)=Transborder_reset(newX(i,:),ub,lb,dim,best_position);
        tmp_fitness=fobj(newX(i,:));
        FEs=FEs+1;
        if tmp_fitness<Fitness(i)
            X(i,:)= newX(i,:);
            Fitness(i)=tmp_fitness;
        end
    end
    [fit_min,minfit_idx]=min(Fitness);
    if fit_min<bestFitness
        best_position=X(minfit_idx,:);
        bestFitness=fit_min;
    end
%     Convergence_curve(it)=bestFitness;
%     best_pos=best_position;
    it=it+1;
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