function [BestSol,Convergence_curve]=SDE( SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
% 参数向量 parameters [n,N_iteration,beta_min,beta_max,pCR]
% n为种群规模，N_iteration为迭代次数
% beta_min 缩放因子下界 Lower Bound of Scaling Factor
% beta_max=0.8; % 缩放因子上界 Upper Bound of Scaling Factor
% pCR 交叉概率 Crossover Probability
% 要求输入数据为列向量（矩阵）
para=[SearchAgents_no,MaxFEs,0.2,0.8,0.2];
%% 差分进化（DE）算法
nPop=para(1); % 种群规模 Population Size
MaxIt=para(2); % 最大迭代次数Maximum Number of Iterations
nVar=dim; % 自变量维数，此例需要优化两个参数c和g Number of Decision Variables
VarSize=[1,dim]; % 决策变量矩阵大小 Decision Variables Matrix Size
beta_min=para(3); % 缩放因子下界 Lower Bound of Scaling Factor
beta_max=para(4); % 缩放因子上界 Upper Bound of Scaling Factor
pCR=para(5); %  交叉概率 Crossover Probability
lb=ones(1,dim).*lb; % 参数取值下界
ub=ones(1,dim).*ub; % 参数取值上界
%% 初始化 Initialization
FEs=0;
empty_individual.Position=[]; % 种群初始化
empty_individual.Cost=[]; % 种群目标函数值初始化
BestSol.Cost=inf; % 最优值初始化
pop=repmat(empty_individual,nPop,1); % 将保存种群信息的结构体扩展为结构体矩阵，行数等于种群大小
for i=1:nPop % 遍历每个个体
    pop(i).Position=init_individual(lb,ub,dim,1);% 随机初始化个体    
    pop(i).Cost=fobj(pop(i).Position) ;% 计算个体目标函数值
    FEs=FEs+1;
    if pop(i).Cost<BestSol.Cost % 如果个体目标函数值优于当前最优值
        BestSol=pop(i); % 更新最优值
    end    
end
BestCost=zeros(MaxIt,1); % 初始化迭代最优值
Convergence_curve=[];
it=1;
z1=0.03; % parameter
%% 主循环 DE Main Loop
while FEs<MaxIt
    for i=1:nPop % 遍历每个个体
        x=pop(i).Position; % 提取个体位置
        % 随机选择三个个体以备变异使用
        A=randperm(nPop); % 个体顺序重新随机排列
        A(A==i)=[]; % 当前个体所排位置腾空（产生变异中间体时当前个体不参与）
        a=A(1);
        b=A(2);
        c=A(3);
        % 变异操作 Mutation
        beta=unifrnd(beta_min,beta_max,VarSize); % 随机产生缩放因子
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position); % 产生中间体
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
        FEs=FEs+1;
        if NewSol.Cost<pop(i).Cost % 如果新个体优于当前个体
            pop(i)=NewSol; % 更新当前个体
            if pop(i).Cost<BestSol.Cost % 如果当前个体（更新后的）优于最优个体
               BestSol=pop(i); % 更新最优个体
            end
        end
    end 
    
    %% SMA
    for po=1:nPop
        X(po,:)=pop(po).Position;%nest暂时保存position
        AllFitness(po)=pop(po).Cost;
    end
%     [X]=SMA(nPop,X,AllFitness,MaxFEs,FEs,lb,ub,dim,BestSol.Position,BestSol.Cost,z1);
    [X]=Fun_SMA(nPop,X,AllFitness,MaxFEs,FEs,lb,ub,dim,BestSol.Position,BestSol.Cost,z1);
    
    for po=1:nPop
        Flag4ub=X(po,:)>ub;        %返回大于ub的逻辑值
		Flag4lb=X(po,:)<lb;
		X(po,:)=(X(po,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
        temp_fit=fobj(X(po,:));
        FEs=FEs+1;
        if temp_fit<pop(po).Cost
            pop(po).Position=X(po,:);%nest暂时保存position
            pop(po).Cost=temp_fit;
        end
        if temp_fit<BestSol.Cost % 如果当前个体（更新后的）优于最优个体
%             BestSol=pop(i); % 更新最优个体
            BestSol.Position=X(po,:);%nest暂时保存position
            BestSol.Cost=temp_fit;
        end
    end

    
    
    
    
    
    % 保存当前迭代最优个体函数值 Update Best Cost
    BestCost(it)=BestSol.Cost;  
    Convergence_curve(it)=BestSol.Cost;
    it=it+1;
end
Positionbest=BestSol.Position;
bestCVaccuarcy=BestSol.Cost;
end
function x=init_individual(xlb,xub,dim,sizepop)
% 参数初始化函数
% lb：参数下界，行向量
% ub：参数上界，行向量
% dim：参数维度
% sizepop 种群规模
% x：返回sizepop*size(lb,2)的参数矩阵
xRange=repmat((xub-xlb),[sizepop,1]);
xLower=repmat(xlb,[sizepop,1]);
x=rand(sizepop,dim).*xRange+xLower;
end