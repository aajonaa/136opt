%% The original AO.

function [bestFitness,Leader_pos,Convergence_curve, Time]=YKAO2(Nagents,NIter,dim,A,trn,vald,TFid,classifierFhd)
    % function [bestfitness,Convergence_curve]=AO0(N,MaxFEs,lb,ub,dim,fobj) % For myself
    %初始化参�?
    lb=0;
    ub=1;
    tic;
    
    it=1;
    Fitnorm=zeros(1,Nagents);
    Convergence_curve=[];
    GV = zeros(Nagents, dim);         % Initial growth vectors
    %% 种群的初始化
    % 初始化一个个�?
    %X=initialization(Nagents,dim,ub,lb);
    X=initialization(Nagents,dim,1,0)>0.5;
    %计算初始种群的�?�应度�??
    for i=1:Nagents
        GV(i, :) = X(i, :) ./ (ub - lb);
        AllFitness(i)=AccSz2(X(i,:)>0,A,trn,vald,classifierFhd);
        %FEs=FEs+1;
    end
    % 为初始种群排序成新种群，找出�?优个体，并记�?
    [fmin,x]=min(AllFitness);
    %
    newX=zeros(Nagents,dim);
    best=X(x,:);
    bestFitness=fmin;
    %% 主循�?
    while it<=NIter
        
        K= 1-((it)^(1/6)/(NIter)^(1/6));
        %E 是一个随迭代递减的数值，它可以成为控制后期局部开发的权重，调节全�?搜索和局部开�?
        E =1*exp(-4*(it/NIter));
        lamda_t = 0.1 + (0.518 * ((1-(it/NIter)^0.5))); %%%%%
    
        for i=1: Nagents
            Fitnorm(i)= (AllFitness(i)-min(AllFitness))/(max(AllFitness)-min(AllFitness));
            for j=1:dim
                %攻击阶段  
                if rand<K %(Comprehensive elimination phase
                    if rand<0.5
                        newX(i,j) = X(i,j)+E.*X(i,j)*(-1)^it; %% 1
                    else
                        % while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= x, break, end, end
                        % while true, r2 = round(N * rand + 0.5); if r2 ~= i && r2 ~= x && r2 ~= r1, break, end, end
                        % omega_it = rand; %%%%%
                        % newX(i,j) = best(j) + ((X(r2,j) - X(i,j)) * lamda_t) + ((X(r1,j) - X(i,j)) * omega_it);
                        newX(i,j) = X(i,j)+E.*best(j)*(-1)^it; %% 2
                        % ori_value = rand(1,1);
                        % cauchy_value = tan((ori_value-0.5)*pi);
                        % newX(i,j)=X(i,j)+cauchy_value* (X(i,j)-best(j));
                    end
                else
                    % DOF = exp(-(2*sinh(FEs/Max_FEs)^0.5-rand))^2;
                    % newX(i,j)= best(j) + DOF .* ((ub-lb)*rand+lb);
                    % ori_value = rand(1,1);
                    % cauchy_value = tan((ori_value-0.5)*pi);
                    % newX(i,j)=X(i,j)+cauchy_value* (X(i,j)-best(j));
                    newX(i,j)=X(i,j);
                    % while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= x, break, end, end
                    % while true, r2 = round(N * rand + 0.5); if r2 ~= i && r2 ~= x && r2 ~= r1, break, end, end
                    % omega_it = rand; %%%%%
                    % newX(i,j) = best(j) + ((X(r2,j) - X(i,j)) * lamda_t) + ((X(r1,j) - X(i,j)) * omega_it);
                end

                %% Local clearance changed to this
                % if rand<Fitnorm(i) %(Local clearance phase
                    % A=randperm(N);
                    % beta=(rand/2)+0.1;
                    omega_it=(rand/2)+0.1;
                    % newX(i,j)=X(A(3),j)+beta.*(X(A(1),j)-X(A(2),j)); %% 3
                    while true, r1 = round(Nagents * rand + 0.5); if r1 ~= i && r1 ~= x, break, end, end
                    while true, r2 = round(Nagents * rand + 0.5); if r2 ~= i && r2 ~= x && r2 ~= r1, break, end, end
                    % omega_it = rand; %%%%%
                    newX(i,j) = best(j) + ((X(r2,j) - X(i,j)) * lamda_t) + ((X(r1,j) - X(i,j)) * omega_it);
                % end
            end


            %% Post-consolidation changed to this
            ii = i + 1;
            if i == Nagents
                ii = 1;
            end
            beta_1 = 1 + (rand / 2); % beta value in Algorithm 1 (line 8)
            if  AllFitness(i) < beta_1 * bestFitness
                newX(i, :) = X(i, :) + abs(randn(1, dim)) .* (X(ii, :) - X(i, :)) + randn(1, dim) .* GV(i, :);
            end
            GV(i, :) = GV(i, :) .* ((rand ^ 2) * randn(1, dim));

            %% 边界收束 ：重新给药，并不是所有的二氧青蒿素都会寻找到疟原虫，�?些会因为自身代谢排出人体，需要再次给药�??
            
            % newX(i,:)=Mutation(newX(i,:),X(i,:),best,dim); %% 4 %(Post-consolidation phase
            newX(i,:)=Transborder_reset(newX(i,:),ub,lb,dim,best); 
            % 计算适应度�??
            
            
            tFitness=AccSz2((newX(i,:)>0),A,trn,vald,classifierFhd);
            
            
            
            %FEs=FEs+1;
            %初步更新适应度�??
            if tFitness<AllFitness(i)
                X(i,:)= newX(i,:);
                AllFitness(i)=tFitness;
            end
        end
        [fmin,x]=min(AllFitness);
        if fmin<bestFitness
            best=X(x,:);
            bestFitness=fmin;
        end

        % Sort Population by Cost
        [~, SortOrder] = sort(AllFitness);
        
        % Ensure that SortOrder size doesn't exceed array bounds
        SortOrder = SortOrder(1:Nagents);  % Only keep the first N sorted indices

        X = X(SortOrder, :);
        AllFitness = AllFitness(SortOrder);
        GV = GV(SortOrder, :);

        % 保存每次迭代的最佳函数�??
        Convergence_curve(it)=bestFitness;
        Leader_pos=best>0;
        bestFitness = min(AllFitness);
        it=it+1;
    end
    Time = toc;
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