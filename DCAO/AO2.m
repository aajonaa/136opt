%% The original AO.

function [Leader_pos,Convergence_curve]=AO2(N,MaxFEs,lb,ub,dim,fobj)
    % function [bestfitness,Convergence_curve]=AO0(N,MaxFEs,lb,ub,dim,fobj) % For myself
    %初始化参数
    FEs=0;
    it=1;
    Fitnorm=zeros(1,N);
    Convergence_curve=[];
    %% 种群的初始化
    % 初始化一个个体
    pop=initialization(N,dim,ub,lb);
    %计算初始种群的适应度值
    for i=1:N
        Fitness(i)=fobj(pop(i,:));
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
                        New_pop(i,j) = pop(i,j)+E.*pop(i,j)*(-1)^FEs; %% 1
                    else
                        New_pop(i,j) = pop(i,j)+E.*best(j)*(-1)^FEs; %% 2
                    end
                else
                    New_pop(i,j)=pop(i,j);
                end
            end
            for j = 1:dim
                if rand<Fitnorm(i)
                    A=randperm(N);
                    beta=(rand/2)+0.1;
                    New_pop(i,j)=pop(A(3),j)+beta.*(pop(A(1),j)-pop(A(2),j)); %% 3
                end
            end

            % New_pop(i, :) = mutation2(New_pop(i, :), N, dim, pop, Fitness, FEs, MaxFEs);
            
            New_pop(i,:)=Mutation(New_pop(i,:),pop(i,:),best,dim); %% 4
            % New_pop(i, :) = mutation2(New_pop(i, :), N, dim, pop, Fitness, FEs, MaxFEs);
            %% 边界收束 ：重新给药，并不是所有的二氧青蒿素都会寻找到疟原虫，一些会因为自身代谢排出人体，需要再次给药、
            New_pop(i,:)=Transborder_reset(New_pop(i,:),ub,lb,dim,best); 
            % 计算适应度值
            tFitness=fobj(New_pop(i,:));
            FEs=FEs+1;

            trail_X = pop(i, :);
            trail_X = mutation2(trail_X, N, dim, pop, Fitness, FEs, MaxFEs);
            ttFitness = fobj(trail_X);
            FEs = FEs + 1;

            if tFitness < ttFitness
                %初步更新适应度值
                if tFitness<Fitness(i)
                    pop(i,:)= New_pop(i,:);
                    Fitness(i)=tFitness;
                end
            else
                if ttFitness<Fitness(i)
                    pop(i,:)= trail_X;
                    Fitness(i)=ttFitness;
                end
            end
        end
        [fmin,x]=min(Fitness);
        if fmin<bestfitness
            best=pop(x,:);
            bestfitness=fmin;
        end
        % 保存每次迭代的最佳函数值
        Convergence_curve(it)=bestfitness;
        Leader_pos=best;
        bestfitness = min(Fitness);
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

function  z = mutation2(z, N, dim, input_pop, input_fitness, FEs, MaxFEs)
    lamda_t = 0.1 + (0.518 * ((1-(FEs/MaxFEs)^0.5))); %%%%%
    [~, sorted_index] = sort(input_fitness,1,'ascend');
    sorted_population = input_pop(sorted_index,:);
    XX = sorted_population;
    bestInd = 1;
    ngS = 6;
    r1 = round(N * rand + 0.5);
    r2 = ngS + round((N - ngS) * rand + 0.5);
    eta_qKR = zeros(1,N); %%%%%
    phi_qKR = 0.25 + 0.55 * ((0 + ((1:N)/N)) .^ 0.5); %%%%%
    for i = 1:N
        eta_qKR(i) = (round(rand * phi_qKR(i)) + (rand <= phi_qKR(i)))/2; %%%%%
    end
    % Compute learning ability
    omega_it = rand; %%%%% （Comply with the whole poppulation)
    for d = 1:dim
        % if rand <= eta_qKR(i)
            z(d) = XX(bestInd,d) + ((XX(r2,d) - XX(i,d)) * lamda_t) + ((XX(r1,d) - XX(i,d)) * omega_it); %%% 3
        % end
    end
end