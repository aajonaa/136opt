%%

function Leader_pos = AO_Net(N,MaxFEs,dim,fobj)
%初始化参数
FEs=0;
% it=1;
Fitnorm=inf * ones(1,N);
% Convergence_curve=[];
%% 种群的初始化
% 初始化一个个体
lb = 0;
ub = 1;
pop=initialization(N,dim,ub,lb);
%计算初始种群的适应度值
for i=1:N
    Fitness(i)=fobj(pop(i,:));
    disp(['training learningRate:', num2str(pop(i, :))]);
    disp(['crossEntropy:', num2str(Fitness(i))]);
    FEs=FEs+1;
end
% 为初始种群排序成新种群，找出最优个体，并记录
[fmin,x]=min(Fitness);
%
New_pop=zeros(N,dim);
best=pop(x,:);
bestfitness=fmin;
%% 主循环
while FEs<=MaxFEs
    
    K= 1-((FEs)^(1/6)/(MaxFEs)^(1/6));
    %E 是一个随迭代递减的数值，它可以成为控制后期局部开发的权重，调节全局搜索和局部开发
    E =1*exp(-4*(FEs/MaxFEs));

    for i=1: N
        
        Fitnorm(i)= (Fitness(i)-min(Fitness))/(max(Fitness)-min(Fitness));
        for j=1:dim
            %攻击阶段  
            if rand<K
                if rand<0.5
                    New_pop(i,j) = pop(i,j)+E.*pop(i,j)*(-1)^FEs;
                else
                    New_pop(i,j) = pop(i,j)+E.*best(j)*(-1)^FEs;
                end
            else
                New_pop(i,j)=pop(i,j);
            end
            
            if rand<Fitnorm(i)
                A=randperm(N);
                beta=(rand/2)+0.1;
                New_pop(i,j)=pop(A(3),j)+beta.*(pop(A(1),j)-pop(A(2),j));
                
            end
        end
        
        New_pop(i,:)=Mutation(New_pop(i,:),pop(i,:),best,dim);
        %% 边界收束 ：重新给药，并不是所有的二氧青蒿素都会寻找到疟原虫，一些会因为自身代谢排出人体，需要再次给药、
        New_pop(i,:)=Transborder_reset(New_pop(i,:),ub,lb,dim,best);
        % 计算适应度值
        tFitness=fobj(New_pop(i,:));
        disp(['training learningRate:', num2str(New_pop(i,:))]);
        disp(['crossEntropy:', num2str(tFitness)]);
        FEs=FEs+1;
        %初步更新适应度值
        if tFitness<Fitness(i)
            pop(i,:)= New_pop(i,:);
            Fitness(i)=tFitness;
        end
    end
    [fmin,x]=min(Fitness);
    if fmin<bestfitness
        best=pop(x,:);
        bestfitness=fmin;
    end
    % 保存每次迭代的最佳函数值
%     Convergence_curve(it)=bestfitness;
    Leader_pos=best;
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