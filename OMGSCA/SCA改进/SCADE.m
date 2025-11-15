% Nenavath, H. and R. K. Jatoth (2018). "Hybridizing sine cosine algorithm with differential evolution for global optimization and object tracking." Applied Soft Computing Journal 62: 1019-1043.

function [Destination_position, Convergence_curve]=SCADE(N,MaxFEs,lb,ub,dim,fobj)
fes = 0;
tic
%%%%***DE中的参数，按照该作者的设置***%%%
beta_min=0.2; % 缩放因子下界 Lower Bound of Scaling Factor
beta_max=0.8; % 缩放因子上界 Upper Bound of Scaling Factor
pCR=0.8; %  交叉概率 Crossover Probability
VarSize=[1,dim];
%%%%***DE中的参数，按照该作者的设置***%%%
%Initialize the set of random solution

X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Objective_values = zeros(1,size(X,1));

% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
    fes = fes + 1;
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
end

%Main loop
t=1; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while fes<=MaxFEs
    
    % Eq. (3.4)
    a = 2;
    r1=a-fes*((a)/MaxFEs); % r1 decreases linearly from a to 0
    
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            
            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
            
        end
    end
    %%*****Muation-Crossover****%%%
     for i=1:1:size(X,1) % 遍历每个个体
        x=X(i,:); % 提取个体位置
  
        % 随机选择三个个体以备变异使用
        A=randperm(size(X,1)); % 个体顺序重新随机排列
        A(A==i)=[]; % 当前个体所排位置腾空（产生变异中间体时当前个体不参与）
        a=A(1);
        b=A(2);
        c=A(3);
        % 变异操作 Mutation
        beta=unifrnd(beta_min,beta_max,VarSize); % 随机产生缩放因子
        y=X(a)+beta.*(X(b)-X(c)); % 产生中间体
        % 防止中间体越界
        y=max(y,lb);
		y=min(y,ub);
        % 交叉操作 Crossover
        z=zeros(size(x)); % 初始化一个新个体
        j0=randi([1,numel(x)]); % 产生一个伪随机数，即选取待交换维度编号？？？
      
        for j=1:numel(x) % 遍历每个维度
            if j==j0 || rand<=pCR % 如果当前维度是待交换维度或者随机概率小于交叉概率
                z(j)=y(j); % 新个体当前维度值等于中间体对应维度值
            else
                z(j)=x(j); % 新个体当前维度值等于当前个体对应维度值
            end
        end
       
        NewSol.Position=z; % 交叉操作之后得到新个体
        NewSol.Cost=fobj(NewSol.Position); % 新个体目标函数值
        fes = fes + 1;
        if NewSol.Cost<fobj(X(i,:)) % 如果新个体优于当前个体
           X(i,:)=NewSol.Position; % 更新当前个体
%             if pop(i).Cost<BestSol.Cost % 如果当前个体（更新后的）优于最优个体
%                BestSol=pop(i); % 更新最优个体
%             end
        end
    end    
   %%*****Muation-Crossover****%%%
    
    for i=1:size(X,1)
         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate the objective values
        Objective_values(1,i)=fobj(X(i,:));
        fes = fes + 1;
        
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    
    Convergence_curve(t)=Destination_fitness;
    
    % Increase the iteration counter
    t=t+1;
end
toc
end
