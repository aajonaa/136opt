%% The original AO.

function [Leader_pos,Convergence_curve]=DCAO(N,MaxFEs,Lb,Ub,dim,fobj)
    % function [bestfitness,Convergence_curve]=AO0(N,MaxFEs,lb,ub,dim,fobj) % For myself
    %初始化参数
    FEs=0;
    it=1;
    Fitnorm=zeros(1,N);
    Convergence_curve=[];
    %% 种群的初始化
    % 初始化一个个体
    X=initialization(N,dim,Ub,Lb);
    %计算初始种群的适应度值
    for i=1:N
        AllFitness(i)=fobj(X(i,:));
        FEs=FEs+1;
    end
    % 为初始种群排序成新种群，找出最优个体，并记录
    [fmin,x]=min(AllFitness);
    %
    newX=zeros(N,dim);
    best_pos=X(x,:);
    bestfitness=fmin;
    %% 主循环
    while FEs<=MaxFEs
        
        K= 1-((FEs)^(1/6)/(MaxFEs)^(1/6));
        %E 是一个随迭代递减的数值，它可以成为控制后期局部开发的权重，调节全局搜索和局部开发
        E =1*exp(-4*(FEs/MaxFEs));
    
        % Sort population by fitness values
        [X, AllFitness, ~] = PopSort(X,AllFitness);
        % Reset
        bestInd = 1;

        for i=1: N
            
            Fitnorm(i)= (AllFitness(i)-min(AllFitness))/(max(AllFitness)-min(AllFitness));

            r1 = round(N * rand + 0.5);
            r2 = 6 + round((N - 6) * rand + 0.5);

            for j=1:dim
                %攻击阶段  
                if rand<K
                    % if rand<0.5
                    %     newX(i,j) = X(i,j)+E.*X(i,j)*(-1)^FEs; %% 1
                    % else
                    %     newX(i,j) = X(i,j)+E.*best_pos(j)*(-1)^FEs; %% 2
                    % end


                    omega_it = rand();
                    lamda_t = 0.1 + (0.518 * ((1-(FEs/MaxFEs)^0.5)));
                    newX(i,j) = X(bestInd, j) + ((X(r2,j) - X(i,j)) * lamda_t) + ((X(r1,j) - X(i,j)) * omega_it);

                else
                    % newX(i,j)=X(i,j);
                    newX(i,j) = X(r1,j) + LnF3(0.518,0.05,1,1);
                end
                
                % if rand<Fitnorm(i)
                %     % A=randperm(N);
                %     % beta=(rand/2)+0.1;
                %     % newX(i,j)=X(A(3),j)+beta.*(X(A(1),j)-X(A(2),j)); %% 3
                %     newX(i,:) = Lb + rand * (Ub - Lb);
                % end
            end
            
            newX(i,:)=Mutation(newX(i,:),X(i,:),best_pos,dim); %% 4
            %% 边界收束 ：重新给药，并不是所有的二氧青蒿素都会寻找到疟原虫，一些会因为自身代谢排出人体，需要再次给药、
            newX(i,:)=Transborder_reset(newX(i,:),Ub,Lb,dim,best_pos); 
            % 计算适应度值
            tFitness=fobj(newX(i,:));
            FEs=FEs+1;
            %初步更新适应度值
            if tFitness<AllFitness(i)
                X(i,:)= newX(i,:);
                AllFitness(i)=tFitness;
            end
        end
        [fmin,x]=min(AllFitness);
        if fmin<bestfitness
            best_pos=X(x,:);
            bestfitness=fmin;
        end
        % 保存每次迭代的最佳函数值
        Convergence_curve(it)=bestfitness;
        Leader_pos=best_pos;
        bestfitness = min(AllFitness);
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

function [sorted_population, sorted_fitness, sorted_index] = PopSort(input_pop,input_fitness)
    [sorted_fitness, sorted_index] = sort(input_fitness,1,'ascend');
    sorted_population = input_pop(sorted_index,:);
end

function Y = LnF3(alpha, sigma, m, n)
    Z = laplacernd(m, n);
    Z = sign(rand(m,n)-0.5) .* Z;
    U = rand(m, n);
    R = sin(0.5*pi*alpha) .* tan(0.5*pi*(1-alpha*U)) - cos(0.5*pi*alpha);
    Y = sigma * Z .* (R) .^ (1/alpha);
end

function x = laplacernd(m, n)
    u1 = rand(m, n);
    u2 = rand(m, n);
    x = log(u1./u2);
end