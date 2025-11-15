function [X,fes,count,conv] = CPA_func_C(N,Max_iteration,X,Target,Target_fitness,lb,ub,dim,fes,fobj)
count=0;
conv=[];

     a = exp(9-18*fes/Max_iteration);%公式5更新a   
    S0=a*(1-fes/Max_iteration); % r1 decreases linearly from a to 0

     %%%%CMAES Parameters
            XX = X;
            m = dim;
            n = N;      
            mean_X = Target';
            Sum = sum(XX);
            x_ba = (1/n) .*Sum;          
            for i=1:n
               dis(i) = sum((XX(i,:) - x_ba).^(2));
            end            
            d = sum((XX(1,:) - x_ba).^2);
            cita = (d)/(sum(dis));          
             mu = n/2;%%%%%%%%%%%%%%%%%%%%%%%%%%%注意
            %sigma = 0.05;
            weights = log(mu+1/2)-log(1:mu)'; % muXone recombination weights muXone复合权重
            mu = floor(mu); % number of parents/points for recombination
            weights = weights/sum(weights); % normalize recombination weights array 规格化重组权值数绿
            mueff=sum(weights)^2/sum(weights.^2);
%             size(weights)
            cc = (4+mueff/m) / (m+4 + 2*mueff/m); % time constant for cumulation for C C的累积时间常敿 进化路径Pc的学习率
            cs = (mueff+2)/(m+mueff+5); % t-const for cumulation for sigma control t-const用于西格玛控制的累积  cs为步长的学习玿
            c1 = 2 / ((m+1.3)^2+mueff); % learning rate for rank-one update of C C皿        rank-1-update的更新学习率
            cmu = 2 * (mueff-2+1/mueff) / ((m+2)^2+2*mueff/2); % and for                     rank-mu-update 的更新学习率
            damps = 1 + 2*max(0, sqrt((mueff-1)/(m+1))-1) + cs;
            pc = zeros(m,1); %进化路径
            ps = zeros(m,1); % evolution paths for C and sigma C和sigma的演化路徿  ps为步长进化路徿
            B = eye(m); % B defines the coordinate system B定义了坐标系  B为一个N*N的单位阵
            D = eye(m); % diagonal matrix D defines the scaling 对角矩阵D定义了缩政 D为一个N*N的单位阵（注意变量用重复）
            C = B*D*(B*D)'; % covariance matrix  C为协方差矩阵
            chiN=m^0.5*(1-1/(4*m)+1/(21*m^2));
            %%%%%%%%
    
    
    
    
    
    
    %Communication and Collaboration
    for j=1 : dim
        EG=X(:,j);     %Selecting a dimensional population
        EG=sort(EG);
        i=round(rand()*(size(X,1))+0.5);  %Random selection of an independent individual
        if rand()< fes/Max_iteration
             X(i,:)=Target;
        end
            r=rand();
            X(i,j)=r*X(i,j)+(1-r)*(EG(1)+EG(2))/2; %Improve its position with two best individuals
    end
    

    for i=1:N % in i-th solution
        

        
        S=2*S0*rand-S0; %A represents the strength of the prey, which decreases with the number of iterations
        
        l=rand();
        
        if  abs(S) < 2*a/3 
            if  rand() > 0.5
                %Disperse food
                X(i,:)=Target-S*(rand(1,dim)*(ub(1)-lb(1))+lb(1));
            else
                %Encircle food
                Dp=abs(Target-X(i,:));
                X(i,:)=Target-2*S*Dp*exp(l)*tan(l*pi/4);
            end
        else
            
            distance=4*rand()-2;   %Represents distance from prey
            
            if abs(distance) < 1  %% [-2) -1...1 (2]  公式14
                %Supporting closest individual
                for j= 1:dim
                    P(j,:)=rand(1,dim)*X(i,j);
                    p_value(j,:)=fobj(P(j,:));
                    fes=fes+1;
                    count = count + 1;
                    conv(1,count) = Target_fitness;
                end
                [p_value,index]=sort(p_value);
                X(i,:)=P(index(1),:);
            else                                       %公式14
                %Research food
                X_rand=rand(1,dim).*(ub(1)-lb(1))+lb(1);
                DD=abs(2*rand()*X_rand-X(i,:));
                %%
                X(i,:)=X_rand-S*DD;

            end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%****cma***%%%%
    [CMA_X,mean_X,ps,pc,C,cita,B,D,fes] = CMS(mean_X,cita,mu,fes,cc,cs,c1,cmu,damps,pc,ps,B,D,C,chiN,fobj,n,dim,weights,mueff);
    CMA_X =CMA_X';
for ii=1:size(CMA_X,1)
    fes = fes + 1;
    Temp_fit(ii) = fobj(CMA_X(ii,:));
    if Temp_fit(ii) < Target_fitness
       Target_fitness = Temp_fit(ii);
       Target = CMA_X(ii,:);
    end
end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        

    end
    
    


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

