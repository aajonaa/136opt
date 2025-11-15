

% CMAES
function [best_score,Convergence_curve]=CMAESS(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)%
% 此处N指问题维度  lambda是指种群大小。
% initialize alpha, beta, and delta_pos


xmean=zeros(1,dim);
best_score=inf; %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Positions = border_control(Positions,lb,ub);
for i=1:SearchAgents_no
    fitness_Val(i) = fobj(Positions(i,:));
end
[~,sort_index] = sort(fitness_Val);
xmean = Positions(sort_index(1),:)';
display(xmean);
best_score = fitness_Val(sort_index(1));
noise = 0.1;
N = dim; % number of objective variables/problem dimension 目标变量数量/问题维度
sigma = 0.05; % coordinate wise standard deviation (step-size) 坐标标准差(步长)
stopeval = MaxFEs;

% Strategy parameter setting: Selection 策略参数设置:选择
%lambda = 4+floor(3*log(dim)); % population size, offspring number 种群大小，后代数量 种群标准差更新时的选中λ个个体
lambda = SearchAgents_no;
mu = lambda/2; % lambda=12; mu=3; weights = ones(mu,1); would be (3_I,12)-ES
weights = log(mu+1/2)-log(1:mu)'; % muXone recombination weights muXone复合权重
mu = floor(mu); % number of parents/points for recombination
weights = weights/sum(weights); % normalize recombination weights array 规格化重组权值数组
mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu 变量-有效的mu大小   方差影响选择集

% Strategy parameter setting: Adaptation 策略参数设置:自适应
cc = (4+mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C C的累积时间常数  进化路径Pc的学习率
cs = (mueff+2)/(N+mueff+5); % t-const for cumulation for sigma control t-const用于西格玛控制的累积  cs为步长的学习率
c1 = 2 / ((N+1.3)^2+mueff); % learning rate for rank-one update of C C的         rank-1-update的更新学习率
cmu = 2 * (mueff-2+1/mueff) / ((N+2)^2+2*mueff/2); % and for                     rank-mu-update 的更新学习率
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma          阻尼系数为dσ

% Initialize dynamic (internal) strategy parameters and constants 初始化动态(内部)策略参数和常量
pc = zeros(N,1); %进化路径
ps = zeros(N,1); % evolution paths for C and sigma C和sigma的演化路径   ps为步长进化路径
B = eye(N); % B defines the coordinate system B定义了坐标系  B为一个N*N的单位阵
D = eye(N); % diagonal matrix D defines the scaling 对角矩阵D定义了缩放  D为一个N*N的单位阵
C = B*D*(B*D)'; % covariance matrix  C为协方差矩阵
%eigeneval = 0; % B and D updated at counteval == 0 B和D在counteval == 0时更新
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2)); % expectation of  ……的期望    chiN 指E||N(0,I)||

%Opt_Record = [];
Convergence_curve=[];
n=0;
counteval = 1;
FEs = 0;
% Main loop
while FEs < stopeval

   % Generate and evaluate lambda offspring 生成和计算lambda后代
    for k=1:lambda
        arz(:,k) = randn(N,1); % standard normally distributed vector 标准正态分布向量
        arx(:,k) = xmean + sigma * (B*D * arz(:,k)); % add mutation % Eq. 40 添加突变 采样操作
		FEs = FEs + 1;
        arfitness(k) = fobj(arx(:,k)') + noise*randn; % objective function call 目标函数调用
    end
    
    %% 绘制响应面和采集函数
    %subplot(1,2,2);
    %mesh(x_o,x_o,y_o);
    %hold on
    %plot3(arx(2,:),arx(1,:),arfitness,'k.','MarkerSize',20);
    %plot3(xmean(2),xmean(1),target_function(xmean),'r.','MarkerSize',30);
    %axis([0 3 0 3 0 9])
    %title(['Iteration = ' num2str(counteval+1)])
    %hold off
    n = n+1;
    %saveas(gcf,['baysian_' num2str(n) '.bmp'])

    % Sort by fitness and compute weighted mean into xmean 排序的适应度和计算加权平均成xmean
    [arfitnesss, arindex] = sort(arfitness); % minimization
    temp_score = arfitness(arindex(1));
    if temp_score < best_score
        best_score = temp_score;
    end
    xmean = arx(:,arindex(1:mu))*weights; % recombination % Eq. 42  更新种群期望值
    zmean = arz(:,arindex(1:mu))*weights; % == D?-1*B’*(xmean-xold)/sigma

    % Cumulation: Update evolution paths 累积:更新进化路径
    ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B * zmean); % Eq. 43  
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*FEs/lambda))/chiN < 1.4+2/(N+1); % 这个具体指什么不太懂 chiN 指E||N(0,I)||
    pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (B*D*zmean); % Eq. 45

    % Adapt covariance matrix C 调整协方差矩阵C
    C = (1-c1-cmu) * C ... % regard old matrix % Eq. 47
        + c1 * (pc*pc' ... % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction
        + cmu ... % plus rank mu update
        * (B*D*arz(:,arindex(1:mu))) ...
        * diag(weights) * (B*D*arz(:,arindex(1:mu)))';

    % Adapt step-size sigma 调整步长σ
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); % Eq. 44

    % Update B and D from C 从C更新B和D
    C=triu(C)+triu(C,1)'; % enforce symmetry 执行对称
    [B,D] = eig(C); % eigen decomposition, B==normalized eigenvectors 特征分解，B==归一化特征向量 B 是协方差矩阵C的特征值构成的矩阵，D是协方差矩阵C的特征向量构成的矩阵
    D = diag(sqrt(diag(D))); % D contains standard deviations now D包含标准差

    %Opt_Record = [Opt_Record xmean];
    Convergence_curve(counteval)=best_score;
    counteval = counteval + 1;
 

end % while, end generation loop
end

function V = border_control(X,L,U)
for i = 1:size(X,1)
    Flag4ub=X(i,:)>U;
    Flag4lb=X(i,:)<L;
    V(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+U.*Flag4ub+L.*Flag4lb;
end
end




