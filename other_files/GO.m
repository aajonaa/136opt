% Jona Mod 2025-1-12
% Growth Optimizer: A powerful metaheuristic algorithm for solving different optimization problems
function [best_pos,Convergence_curve]= GO(N, MaxFEs, lb, ub, dim, fobj)
FEs=0;
X=lb+(ub-lb)*unifrnd(0,1,N,dim);
Convergence_curve = [];
it = 1;
bestFitness=inf;
for i=1:N
    AllFitness(i)=fobj(X(i,:));
    FEs=FEs+1;
    if bestFitness>AllFitness(i)
        bestFitness=AllFitness(i);
        best_pos=X(i,:);
    end
end
Convergence_curve(it)=bestFitness;
while FEs <= MaxFEs
    [~, ind]=sort(AllFitness);
    Best_X=X(ind(1),:);
%% Learning phase
    for i=1:N
        Worst_X = X(ind(randi([N-4,N])),:);
        Better_X=X(ind(randi([2,5])),:);
        random=selectID(N,i,2);
        L1=random(1);
        L2=random(2);
        D_value1=(Best_X-Better_X);
        D_value2=(Best_X-Worst_X);
        D_value3=(Better_X-Worst_X);
        D_value4=(X(L1,:)-X(L2,:));
        Distance1=norm(D_value1);
        Distance2=norm(D_value2);
        Distance3=norm(D_value3);
        Distance4=norm(D_value4);
        rate=Distance1+Distance2+Distance3+Distance4;
        LF1=Distance1/rate;
        LF2=Distance2/rate;
        LF3=Distance3/rate;
        LF4=Distance4/rate;
        SF=(AllFitness(i)/max(AllFitness));
        Gap1=LF1*SF*D_value1;
        Gap2=LF2*SF*D_value2;
        Gap3=LF3*SF*D_value3;
        Gap4=LF4*SF*D_value4;
        newX(i,:)=X(i,:)+Gap1+Gap2+Gap3+Gap4;
        %Clipping
        newX(i,:)=max(newX(i,:),lb);
        newX(i,:)=min(newX(i,:),ub);
        newAllFitness=fobj(newX(i,:));
        FEs=FEs+1;
        %Update
        if AllFitness(i)>newAllFitness
            AllFitness(i)=newAllFitness;
            X(i,:)=newX(i,:);
        else
            if rand<0.001&&ind(i)~=1
                AllFitness(i)=newAllFitness;
                X(i,:)=newX(i,:);
            end
        end
        if bestFitness>AllFitness(i)
            bestFitness=AllFitness(i);
            best_pos=X(i,:);
        end
    end

%% Reflection phase
    for i=1:N
        newX(i,:)=X(i,:);
        j=1;
        while j<=dim
            if rand<0.3
                R=X(ind(randi(5)),:);
                newX(i,j) = X(i,j)+(R(:,j)-X(i,j))*unifrnd(0,1);
                if rand<(0.01+(0.1-0.01)*(1-FEs/MaxFEs))
                    newX(i,j)=lb+(ub-lb)*unifrnd(0,1);
                end
            end
            j=j+1;
        end
        %Clipping
        newX(i,:)=max(newX(i,:),lb);
        newX(i,:)=min(newX(i,:),ub);
        newAllFitness=fobj(newX(i,:));
        FEs=FEs+1;
        %Update
        if AllFitness(i)>newAllFitness
            AllFitness(i)=newAllFitness;
            X(i,:)=newX(i,:);
        else
            if rand<0.001&&ind(i)~=1
                AllFitness(i)=newAllFitness;
                X(i,:)=newX(i,:);
            end
        end
        if bestFitness>AllFitness(i)
            bestFitness=AllFitness(i);
            best_pos=X(i,:);
        end
    end
    it = it + 1;
    Convergence_curve(it) = bestFitness;
end
end

%Subfunction:  Select four individuals different from Xi
function random = selectID(N, i, Num)
    if Num <= N
        vec=[1:i-1,i+1:N];
        r=zeros(1,2);
        for k =1:2
            n = N-k;
            t = randi(n,1,1);
            r(k) = vec(t);
            vec(t)=[];
        end
        random(1)=r(1);random(2)=r(2);
    end
end

