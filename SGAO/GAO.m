% SGAO (Social-Growth Aquila Optimizer) is an enhanced Aquila algorithm that integrates social learning strategies for knowledge sharing and growth-pattern mechanisms for adaptive exploration, creating a balanced optimization approach.

% Reference1: Poomin Duankhan, Khamron Sunat, Sirapat Chiewchanwattana, Patchara Nasa-ngium,
% The Differentiated Creative Search (DCS): Leveraging differentiated knowledge-acquisition and creative realism to address complex optimization problems
% KBS journal

% Reference2: Mojtaba Ghasemi, Mohsen Zare, Pavel Trojovský, Ravipudi Venkata Rao, Eva Trojovská, Venkatachalam Kandasamy,
% Optimization based on the smart behavior of plants with its engineering applications: Ivy algorithm
% KBS journal
function [best_pos,Convergence_curve]=SGAO(N,MaxFEs,lb,ub,dim,fobj)
    %% 初始化参数
    FEs=0;
    it=1;
    Fitnorm=zeros(1,N);
    Convergence_curve=[];
    GV = zeros(N, dim);         % Initial growth vectors
    %% 种群的初始化
    % 初始化一个个体
    X=initialization(N,dim,ub,lb);
    %计算初始种群的适应度值
    for i=1:N
        GV(i, :) = X(i, :) ./ (ub - lb); 
        AllFitness(i)=fobj(X(i,:));
        FEs=FEs+1;
    end
    %% 为初始种群排序成新种群，找出最优个体，并记录
    [fmin,idx]=min(AllFitness);

    newX=zeros(N,dim);
    Best_pos=X(idx,:);
    bestFitness=fmin;
    %% 主循环
    while FEs<=MaxFEs
        
        K= 1-((FEs)^(1/6)/(MaxFEs)^(1/6));
        %E 是一个随迭代递减的数值，它可以成为控制后期局部开发的权重，调节全局搜索和局部开发
        E =1*exp(-4*(FEs/MaxFEs));
        lamda_t = 0.1 + (0.518 * ((1-(FEs/MaxFEs)^0.5))); %%%%%
    
        for i=1: N
            Fitnorm(i)= (AllFitness(i)-min(AllFitness))/(max(AllFitness)-min(AllFitness));
            for j=1:dim
                %攻击阶段  
                if rand<K %(Comprehensive elimination phase
                    if rand<0.5
                        newX(i,j) = X(i,j)+E.*X(i,j)*(-1)^FEs; %% 1
                    else
                        newX(i,j) = X(i,j)+E.*Best_pos(j)*(-1)^FEs; %% 2
                    end
                else
                    newX(i,j)=X(i,j);
                end

                %% Local clearance changed to this
                if rand<Fitnorm(i) %(Local clearance phase
                    A=randperm(N);
                    beta=(rand/2)+0.1;
                    newX(i,j)=X(A(3),j)+beta.*(X(A(1),j)-X(A(2),j)); %% 3
                end
                % omega_it=(rand/2)+0.1;
                % while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= idx, break, end, end
                % while true, r2 = round(N * rand + 0.5); if r2 ~= i && r2 ~= idx && r2 ~= r1, break, end, end
                % newX(i,j) = Best_pos(j) + ((X(r2,j) - X(i,j)) * lamda_t) + ((X(r1,j) - X(i,j)) * omega_it);
            end


            %% Post-consolidation changed to this
            % newX(i,:)=Mutation(newX(i,:),X(i,:),best,dim); %% 4 %(Post-consolidation phase
            ii = i + 1;
            if i == N
                ii = 1;
            end
            beta_1 = 1 + (rand / 2); % beta value in Algorithm 1 (line 8)
            if  AllFitness(i) < beta_1 * bestFitness
                newX(i, :) = X(i, :) + abs(randn(1, dim)) .* (X(ii, :) - X(i, :)) + randn(1, dim) .* GV(i, :);
            end
            GV(i, :) = GV(i, :) .* ((rand ^ 2) * randn(1, dim));
            
            %% 边界收束 ：重新给药，并不是所有的二氧青蒿素都会寻找到疟原虫，一些会因为自身代谢排出人体，需要再次给药、
            newX(i,:)=Transborder_reset(newX(i,:),ub,lb,dim,Best_pos); 
            % 计算适应度值
            tFitness=fobj(newX(i,:));
            FEs=FEs+1;
            %初步更新适应度值
            if tFitness<AllFitness(i)
                X(i,:)= newX(i,:);
                AllFitness(i)=tFitness;
            end
        end
        [fmin,idx]=min(AllFitness);
        if fmin<bestFitness
            Best_pos=X(idx,:);
            bestFitness=fmin;
        end

        % Sort Population by Cost
        [~, SortOrder] = sort(AllFitness);
        
        % Ensure that SortOrder size doesn't exceed array bounds
        SortOrder = SortOrder(1:N);  % Only keep the first N sorted indices

        X = X(SortOrder, :);
        AllFitness = AllFitness(SortOrder);
        GV = GV(SortOrder, :);

        % 保存每次迭代的最佳函数值
        Convergence_curve(it)=bestFitness;
        best_pos=Best_pos;
        bestFitness = min(AllFitness);
        it=it+1;
    end
end

function z=Transborder_reset(z,ub,lb,dim,best)
    for j=1:dim
        if z(j)>ub || z(j)<lb
            
            z(j)=best(j);
            
        end
    end
end