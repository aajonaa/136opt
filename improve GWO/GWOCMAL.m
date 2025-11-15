
%用beta来初始化种群，前期在levy和GWO中贪心选择最优，后期用CMA
% Grey Wolf Optimizer
function [Alpha_pos,Convergence_curve]=GWOCMAL(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
N = SearchAgents_no;
it=1;FEs=0;
Convergence_curve=[];
R = 0.5;
flag = 1;
for i=1:SearchAgents_no
    Positions(i,:) = Positions(i,:) *(1+betarnd(1.2,1.2));
end
%cov = covariance(Positions);
% Main loop
while  FEs < MaxFEs
    % % Main loop
    % while l<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        tag = 1;
        if FEs<MaxFEs
            FEs=FEs+1;
            % Calculate objective function for each search agent
            fitness=fobj(Positions(i,:));
            
            % Update Alpha, Beta, and Delta
            if fitness<Alpha_score
                Alpha_score=fitness; % Update alpha
                Alpha_pos=Positions(i,:);
                alpha_no = i;
                if tag == 1
                    beta_no1 = alpha_no;
                    delta_no = alpha_no;
                    tag = 0;
                end
            end
            
            if fitness>Alpha_score && fitness<Beta_score
                Beta_score=fitness; % Update beta
                Beta_pos=Positions(i,:);
                beta_no1 = i;
            end
            
            if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score
                Delta_score=fitness; % Update delta
                Delta_pos=Positions(i,:);
                delta_no = i;
            end
        else
            break;
        end
    end
    % Update the Position of search agents including omegas
    if FEs < MaxFEs/2               
        a=2-FEs*((2)/MaxFEs); % a decreases linearly fron 2 to 0
        for i=1:size(Positions,1)        
            for j=1:size(Positions,2)
                
                r1=rand(); % r1 is a random number in [0,1]
                r2=rand(); % r2 is a random number in [0,1]
                
                A1=2*a*r1-a; % Equation (3.3)
                C1=2*r2; % Equation (3.4)
                
                D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
                X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                
                r1=rand();
                r2=rand();
                
                A2=2*a*r1-a; % Equation (3.3)
                C2=2*r2; % Equation (3.4)
                
                D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
                X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2
                
                r1=rand();
                r2=rand();
                
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)
                
                D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
                X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3
%                 w1 = Alpha_score/(Alpha_score + Beta_score + Delta_score);
%                 w2 = Beta_score/(Alpha_score + Beta_score + Delta_score);
%                 w3 = Delta_score/(Alpha_score + Beta_score + Delta_score);
%                 Positions(i,j) = X1 * w1 + X2 * w2 + w3 * X3; 
%                   temp(1,j) =(X1*w1+X2*w2+X3*w3)/3;% Equation (3.7)
%                 temp(1,j)=(X1+X2+X3)/3;% Equation (3.7) 
              Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)           
        end
%            if fobj(temp) <= fobj(temp1)
%                Positions(i,:) = temp;
%            else
%                Positions(i,:) = temp1;
%            end
        end
		X3=Lenvy(Positions,1);
        [Positions,FEs]=Best_Position(Positions,X3,fobj,FEs); 
    else
        X = Positions;
        m = dim;
        n = SearchAgents_no;
        if flag == 1
            mean_alpha = Alpha_pos';
            mean_beta  = Beta_pos';
            mean_delta = Delta_pos';
            
            Sum = sum(X);
            x_ba = (1/n) .*Sum;
            for i=1:n
                dis(i) = sum((X(i,:) - x_ba).^(2));
            end
            d_alpha = sum((X(alpha_no,:) - x_ba).^2);
            cita_alpha = (d_alpha)/(sum(dis));
            d_beta = sum((X(beta_no1,:) - x_ba).^2);
            cita_beta = (d_beta)/(sum(dis));
            d_delta = sum((X(delta_no,:) - x_ba).^2);
            cita_delta = (d_delta)/(sum(dis));
            
            mu = n/2;
            %sigma = 0.05;
            weights = log(mu+1/2)-log(1:mu)'; % muXone recombination weights muXone复合权重
            mu = floor(mu); % number of parents/points for recombination
            weights = weights/sum(weights); % normalize recombination weights array 规格化重组权值数绿
            mueff=sum(weights)^2/sum(weights.^2);
            
            cc = (4+mueff/m) / (m+4 + 2*mueff/m); % time constant for cumulation for C C的累积时间常敿 进化路径Pc的学习率
            cs = (mueff+2)/(m+mueff+5); % t-const for cumulation for sigma control t-const用于西格玛控制的累积  cs为步长的学习玿
            c1 = 2 / ((m+1.3)^2+mueff); % learning rate for rank-one update of C C皿        rank-1-update的更新学习率
            cmu = 2 * (mueff-2+1/mueff) / ((m+2)^2+2*mueff/2); % and for                     rank-mu-update 的更新学习率
            damps = 1 + 2*max(0, sqrt((mueff-1)/(m+1))-1) + cs;
            pc = zeros(m,1); %进化路径
            ps = zeros(m,1); % evolution paths for C and sigma C和sigma的演化路徿  ps为步长进化路徿
            B = eye(m); % B defines the coordinate system B定义了坐标系  B为一个N*N的单位阵
            D = eye(m); % diagonal matrix D defines the scaling 对角矩阵D定义了缩政 D为一个N*N的单位阵
            C = B*D*(B*D)'; % covariance matrix  C为协方差矩阵
            chiN=m^0.5*(1-1/(4*m)+1/(21*m^2));
            %             C_alpha = C;%cov(X(alpha_no,:));
            %             C_beta = C;%cov(X(beta_no,:));
            %             C_delta = C;%cov(X(delta_no,:));
            flag = 0;
            
        end
        [X_alpha,mean_alpha,ps,pc,C,cita_alpha,B,D,FEs] = CMS(mean_alpha,cita_alpha,mu,FEs,cc,cs,c1,cmu,damps,pc,ps,B,D,C,chiN,fobj,n,dim,weights,mueff);
        [X_beta,mean_beta,ps,pc,C,cita_beta,B,D,FEs] = CMS(mean_beta,cita_beta,mu,FEs,cc,cs,c1,cmu,damps,pc,ps,B,D,C,chiN,fobj,n,dim,weights,mueff);
        [X_delta,mean_delta,ps,pc,C,cita_delta,B,D,FEs] = CMS(mean_delta,cita_delta,mu,FEs,cc,cs,c1,cmu,damps,pc,ps,B,D,C,chiN,fobj,n,dim,weights,mueff);
        Temp_pos = [X_alpha';X_beta';X_delta'];
        for i=1:size(Temp_pos,1)
            FEs = FEs + 1;
            Temp_fit(i) = fobj(Temp_pos(i,:));
        end
        [sorted_fit,index]=sort(Temp_fit);
        sorted_pos = Temp_pos(index,:);
        Positions = sorted_pos(1:size(X_delta,1),:);
        temp_score = sorted_fit(1);
        if temp_score < Alpha_score
            Alpha_score = temp_score;
            Alpha_pos = Positions(1,:);
        end
        
        %     Alpha_score = sorted_fit(1);
        %     Alpha_pos = Positions(1,:);
    end
    %     FEs = FEs + 1;
    %     Convergence_curve(FEs)=Alpha_score;
    Convergence_curve(it)=Alpha_score;
    it=it+1;
end
end
function o=L(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end
function [X,xmean,ps,pc,C,sigma,B,D,counteval]= CMS(xmean,sigma,mu,counteval,cc,cs,c1,cmu,damps,pc,ps,B,D,C,chiN,fobj,lambda,N,weights,mueff)
% Generate and evaluate lambda offspring 生成和计算lambda后代
for k=1:lambda
    arz(:,k) = randn(N,1); % standard normally distributed vector 标准正态分布向量
    arx(:,k) = xmean + sigma * (B*D * arz(:,k)); % add mutation % Eq. 40 添加突变 采样操作
    counteval = counteval + 1;
    arfitness(k) = fobj(arx(:,k)'); %+ noise*randn objective function call 目标函数调用   与原来有出入
end
% Sort by fitness and compute weighted mean into xmean 排序的适应度和计算加权平均成xmean
[arfitnesss, arindex] = sort(arfitness); % minimization
%    temp_score = arfitness(arindex(1));
%    if temp_score < best_score
%        best_score = temp_score;
%    end
%    best_score = arfitness(arindex(1));
xmean = arx(:,arindex(1:mu))*weights; % recombination % Eq. 42  更新种群期望值
zmean = arz(:,arindex(1:mu))*weights; % == D?-1*B’*(xmean-xold)/sigma

% Cumulation: Update evolution paths 累积:更新进化路径
ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B * zmean); % Eq. 43
hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4+2/(N+1); % 这个具体指什么不太懂 chiN 指E||N(0,I)||
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

X = arx;
end

function [X,FEs]=Best_Position(X,X3,fobj,FEs)
n=size(X,1);
for i=1:n
	FEs = FEs + 2;
    y_X=fobj(X(i,:)); 
    y_X3=fobj(X3(i,:));
    y_col=[y_X,y_X3];
    [~,k]=min(y_col);
    switch k
        case 1
            X(i,:)=X(i,:);
        case 2
            X(i,:)=X3(i,:);       
    end
end               
end




function X=Lenvy(X0,k)
n=size(X0,1);
for i=1:n
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn*sigma;
    v=randn;
    step=u./abs(v).^(1/beta);
    o=step;
    X1(i,:)=X0(i,:)*(1+k*o);
end
X=X1;
end