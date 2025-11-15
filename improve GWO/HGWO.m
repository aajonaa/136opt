function [parent_Alpha_Position,Convergence_curve]=HGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
% input_train,output_train,input_test,output_test,para=[30,500,0.2,0.8,0.2];]
% 参数向量 parameters [n,N_iteration,beta_min,beta_max,pCR]
% n为种群规模，N_iteration为迭代次数
% beta_min 缩放因子下界 Lower Bound of Scaling Factor
% beta_max=0.8; % 缩放因子上界 Upper Bound of Scaling Factor
% pCR 交叉概率 Crossover Probability
% 要求输入数据为列向量（矩阵）

%% 利用差分进化-灰狼优化混合算法（DE_GWO）选择最佳的SVR参数
nPop=SearchAgents_no; % 种群规模 Population Size
MaxIt=Max_iter; % 最大迭代次数Maximum Number of Iterations
nVar=dim; % 自变量维数，此例需要优化两个参数c和g Number of Decision Variables
VarSize=[1,nVar]; % 决策变量矩阵大小 Decision Variables Matrix Size
beta_min=0.2; % 缩放因子下界 Lower Bound of Scaling Factor
beta_max=0.8; % 缩放因子上界 Upper Bound of Scaling Factor
pCR=0.2; %  交叉概率 Crossover Probability
lb=lb.*ones(1,dim);
ub=ub.*ones(1,dim);
%% 初始化
% 父代种群初始化
parent_Position=initialization(SearchAgents_no,dim,ub,lb);
parent_Val=zeros(nPop,1); % 目标函数值
for i=1:nPop % 遍历每个个体
    parent_Val(i)=fobj(parent_Position(i,:)); % 计算个体目标函数值
end
% 突变种群初始化
mutant_Position=initialization(SearchAgents_no,dim,ub,lb); % 随机初始化位置
mutant_Val=zeros(nPop,1); % 目标函数值
for i=1:nPop % 遍历每个个体
    mutant_Val(i)=fobj(mutant_Position(i,:)); % 计算个体目标函数值
end
% 子代种群初始化
child_Position=initialization(SearchAgents_no,dim,ub,lb); % 随机初始化位置
child_Val=zeros(nPop,1); % 目标函数值
for i=1:nPop % 遍历每个个体
    child_Val(i)=fobj(child_Position(i,:)); % 计算个体目标函数值
end
%% 确定父代种群中的Alpha,Beta,Delta狼
[~,sort_index]=sort(parent_Val); % 父代种群目标函数值排序
parent_Alpha_Position=parent_Position(sort_index(1),:); % 确定父代Alpha狼
parent_Alpha_Val=parent_Val(sort_index(1)); % 父代Alpha狼目标函数值
parent_Beta_Position=parent_Position(sort_index(2),:); % 确定父代Beta狼
parent_Delta_Position=parent_Position(sort_index(3),:); % 确定父代Delta狼
%% 迭代开始
it = 1;
Fes = 0;
Convergence_curve=[];
Convergence_curve(1)=parent_Alpha_Val;
while Fes < MaxIt
    a=2-Fes*((2)/MaxIt); % 对每一次迭代，计算相应的a值，a decreases linearly fron 2 to 0
    % 更新父代个体位置
    for par=1:nPop % 遍历父代个体
        for var=1:nVar % 遍历每个维度            
            % Alpha狼Hunting
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]            
            A1=2*a*r1-a; % 计算系数A
            C1=2*r2; % 计算系数C
            D_alpha=abs(C1*parent_Alpha_Position(var)-parent_Position(par,var));
            X1=parent_Alpha_Position(var)-A1*D_alpha;
            % Beta狼Hunting
            r1=rand();
            r2=rand();            
            A2=2*a*r1-a; % 计算系数A
            C2=2*r2; % 计算系数C
            D_beta=abs(C2*parent_Beta_Position(var)-parent_Position(par,var));
            X2=parent_Beta_Position(var)-A2*D_beta;
            % Delta狼Hunting
            r1=rand();
            r2=rand();
            A3=2*a*r1-a; % 计算系数A
            C3=2*r2; % 计算系数C
            D_delta=abs(C3*parent_Delta_Position(var)-parent_Position(par,var));
            X3=parent_Delta_Position(var)-A3*D_delta;
            % 位置更新，防止越界
            X=(X1+X2+X3)/3;
            X=max(X,lb(var));
            X=min(X,ub(var));
            parent_Position(par,var)=X;
        end
        Fes = Fes + 1;
        parent_Val(par)=fobj(parent_Position(par,:)); % 计算个体目标函数值
    end
	
	if Fes < MaxIt
		Fes = Fes + 1;
		% 产生变异（中间体）种群
		for mut=1:nPop
			A=randperm(nPop); % 个体顺序重新随机排列
			A(A==i)=[]; % 当前个体所排位置腾空（产生变异中间体时当前个体不参与）%%好像应该把i换成mut
			a=A(1);
			b=A(2);
			c=A(3);
			beta=unifrnd(beta_min,beta_max,VarSize); % 随机产生缩放因子  beta 是一个和 VarSize 同型，大小服从 beta_min到beta_max之间均匀分布的随机值
			y=parent_Position(a)+beta.*(parent_Position(b)-parent_Position(c)); % 产生中间体
			% 防止中间体越界
			y=max(y,lb);
			y=min(y,ub);
			mutant_Position(mut,:)=y;
		end
		% 产生子代种群，交叉操作 Crossover
		for child=1:nPop
			x=parent_Position(child,:);
			y=mutant_Position(child,:);
			z=zeros(size(x)); % 初始化一个新个体
			j0=randi([1,numel(x)]); % 产生一个伪随机数，即选取待交换维度编号？？？
			for var=1:numel(x) % 遍历每个维度
				if var==j0 || rand<=pCR % 如果当前维度是待交换维度或者随机概率小于交叉概率
					z(var)=y(var); % 新个体当前维度值等于中间体对应维度值
				else
					z(var)=x(var); % 新个体当前维度值等于当前个体对应维度值
				end
			end
			child_Position(child,:)=z; % 交叉操作之后得到新个体
			child_Val(child)=fobj(z); % 新个体目标函数值
		end
		% 父代种群更新
		for par=1:nPop
			if child_Val(par)<parent_Val(par) % 如果子代个体优于父代个体
				parent_Val(par)=child_Val(par); % 更新父代个体
			end
		end
		% 确定父代种群中的Alpha,Beta,Delta狼
		[~,sort_index]=sort(parent_Val); % 父代种群目标函数值排序
		parent_Alpha_Position=parent_Position(sort_index(1),:); % 确定父代Alpha狼
		parent_Alpha_Val=parent_Val(sort_index(1)); % 父代Alpha狼目标函数值
		parent_Beta_Position=parent_Position(sort_index(2),:); % 确定父代Beta狼
		parent_Delta_Position=parent_Position(sort_index(3),:); % 确定父代Delta狼
		if (it>1)&(Convergence_curve(it-1)<parent_Alpha_Val)%防止震荡
			Convergence_curve(it)=Convergence_curve(it-1);
		else 
			Convergence_curve(it)=parent_Alpha_Val;
			it = it + 1;
		end
	else
		break;
	end;
end
bestc=parent_Alpha_Position(1,1);
bestg=parent_Alpha_Position(1,2);

end
