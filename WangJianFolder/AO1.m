%%

function [best_pos,Convergence_curve]=AO(N,Max_FEs,lb,ub,dim,fobj) 
% function [bestFitness,Convergence_curve]=AO(N,Max_FEs,lb,ub,dim,fobj) % For myself

FEs=0;
it=1;
Fitnorm=zeros(1,N);
Convergence_curve=[];
FitnessAll = ones(1, N) * inf;

X=initialization(N,dim,ub,lb);

for i=1:N
    FitnessAll(i)=fobj(X(i,:));
    FEs=FEs+1;
end

[fit_min,minfit_idx]=min(FitnessAll);
%
newX=zeros(N,dim);
best_position=X(minfit_idx,:);
bestFitness=fit_min;

while FEs<=Max_FEs
    
    K= 1-((FEs)^(1/3)/(Max_FEs)^(1/3));
    E =1*exp(-4*(FEs/Max_FEs));
    
    %% Traveling distance rate, the original constant is 1/6
    TDR = K;

    for i=1: N
        
        Fitnorm(i)= (FitnessAll(i)-min(FitnessAll))/(max(FitnessAll)-min(FitnessAll));
        for j=1:dim
 
            if rand<K
                if rand<0.5
                    newX(i,j) = X(i,j)+E * X(i,j)*(-1)^FEs; %% 1
                else
                    newX(i,j) = X(i,j)+E * best_position(j)*(-1)^FEs; %% 2
                end
            else
%                 newX(i,j)=X(i,j); %% 3
                %% First part change: Add the sine-learning mechanism (SLM) based on von Mises distribution (VMD)
                theta = pi/2 * min(Fitnorm)/Fitnorm(i);
                newX(i,j)=X(i,j) + tan(theta) * best_position(j); %% 3
            end
            
            if rand<Fitnorm(i)
                A=randperm(N);
                beta=(rand/2)+0.1;
                newX(i,j)=X(A(3),j)+beta *(X(A(1),j)-X(A(2),j)); %% 4
                
            end
        end
        
        newX(i,:)=Mutation(newX(i,:),X(i,:),best_position,dim); %% 5
        newX(i,:)=Transborder_reset(newX(i,:),ub,lb,dim,best_position); %% Boundary control
        tmp_fitness=fobj(newX(i,:));
        FEs=FEs+1;
        if tmp_fitness<FitnessAll(i)
            X(i,:)= newX(i,:);  %% Greedy selection
            FitnessAll(i)=tmp_fitness;
        end
    end
    [fit_min,minfit_idx]=min(FitnessAll);
    if fit_min<bestFitness
        best_position=X(minfit_idx,:);
        bestFitness=fit_min;
    end
    Convergence_curve(it)=bestFitness;
    best_pos=best_position;
    bestFitness = min(FitnessAll);
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