%___________________________________________________________________%
%  Grey Wold Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% Grey Wolf Optimizer
function [Alpha_pos,Convergence_curve]=SOGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
Positions1=initialization(SearchAgents_no,dim,ub,lb);
%Pos_Temp = initialization(SearchAgents_no,dim,ub,lb);
 %for i=1:size(Positions,1)
%        for j=1:size(Positions,2) 
  %          r = exp(j);
 %           Pos_Temp(i,j) = i * r - floor(i * r);
        %    Positions(i,j) =lb + Pos_Temp(i,j) * (ub-lb);% 
       % end
 %end
 FEs=0;
 K = SearchAgents_no;


% Convergence_curve=zeros(1,Max_iter);
 Convergence_curve=[];
 
%%   混沌初始
% x = rand();
% while(~(x~=0.25 && x~=0.5 && x~=0.75 && x~=1))   %x==0.25|x==0.5|x==0.75|x==1
%      x=rand();
% end
%    ch(1) = x;
% for ant=1:(SearchAgents_no-1)
%     ch(ant+1)=4*ch(ant)*(1-ch(ant));    %式（4-1），ch即为该式中的βi
%     PCh = ch(ant)*Positions1;         %式（4-2）
%     PHe = [Positions1;PCh];%phe为混沌变化后的搜索代理位置和变化前的的位置
%     count=size(PHe,1);%count为所有搜索代理的个数
%     FitnessHe1=[];
%     for i=1:count%将所有适应度值计算存放在 FitnessHe1 中
%          FEs=FEs+1;
%          PHeLin=fobj(PHe(i,:));
%          FitnessHe1 = [FitnessHe1 PHeLin];
%     end
% %        FitnessHe1=calculateFitness(ObjValHe1);
%          [FitnessHe2,index] = sort(FitnessHe1);%FitnessHe2保存的是排序后的适应度值，index是FitnessHe2中保存的FitnessHe1的索引
%           Positions = PHe(index,:);
%          Positions = Positions1(1:SearchAgents_no,:);%将 Positions按从小到大的顺序排序
% end
%%




it=1;% Loop counter
q=2;
% Main loop
while FEs<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        if FEs<Max_iter
            FEs=FEs+1;
            % Calculate objective function for each search agent
            fitness=fobj(Positions(i,:));
            
            % Update Alpha, Beta, and Delta
            if fitness<Alpha_score
                Alpha_score=fitness; % Update alpha
                Alpha_pos=Positions(i,:);
            end
            
            if fitness>Alpha_score && fitness<Beta_score
                Beta_score=fitness; % Update beta
                Beta_pos=Positions(i,:);
            end
            
            if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score
                Delta_score=fitness; % Update delta
                Delta_pos=Positions(i,:);
            end
        else
            break;
        end
    end
    
    
%     a=2-FEs*((2)/Max_iter); % a decreases linearly fron 2 to 0
%     a=2-FEs*((2)/Max_iter);% a decreases linearly fron 2 to 0
    n=1-FEs*((1)/Max_iter); 
    
    if FEs <= Max_iter/2
        a = 1 + (2 - 0)*(1 + (cos(((FEs-1)*pi)/(Max_iter-1)))^n)/2;
    else
        a = 1 + (2 - 0)*(1 - abs(cos(((FEs-1)*pi)/(Max_iter-1)))^n)/2;
    end
    
    %%
         %交叉
         for iter=1:q
             for i=1:size(Positions,1)%遍历整个解空间的位置，找到下一个合适的搜索代理的位置
                 for j=1:size(Positions,2)
                     std = floor((log(FEs)/FEs)) * (Positions(i,:) - Alpha_pos);
                     keci = rand();
                     keci1 = rand();
                     Positions(i,:) =(keci * Alpha_pos - keci1 * Positions(i,:));%Positions(i,:) + 11 gaussian(Alpha_score,std)+
                 end
             end
         end
    %%
%     Update the Position of search agents including omegas
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
            
           % Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
            
            w1 = abs(X1)/abs(X1+X2+X3);
            w2 = abs(X2)/abs(X1+X2+X3);
            w3 = abs(X3)/abs(X1+X2+X3);
            Positions(i,j)=(X1*w1+X2*w2+X3*w3)/3;% Equation (3.7)
%              Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)

        end
        
        if(i<SearchAgents_no)
           %%     变异
            FEs=FEs+2;
            Moth_pos_m_gaus=Positions(i,:)*(1+randn(1));%经过高斯变异后搜索代理的位置
            Moth_fitness_m_gaus=fobj(Moth_pos_m_gaus);%经过高斯变异后搜索代理的适应度值
            Moth_fitness_s=fobj(Positions(i,:));%没有经过高斯变异的搜索代理的适应度值
            
            Moth_fitness_comb=[Moth_fitness_m_gaus,Moth_fitness_s];%将两次的适应度值拼接在一起
            [~,m]=min(Moth_fitness_comb);%找到所有适应度中的最小值，m为最小值的下标
            if m==1%如果m为1则表示经过高斯变异后的位置比较好，用高斯变异后的位置取代原来的位置，否则就原来的位置不变
                Positions(i,:)=Moth_pos_m_gaus;
            else
                Positions(i,:)=Positions(i,:);
            end
  %%          
        end   
    end
    
    %%  %使用搜索策略预形成混沌局部搜索
    setCan = (Max_iter-FEs+1)/Max_iter;%相当于入
    x = rand();
    while(~(x~=0.25 && x~=0.5 && x~=0.75 && x~=1))
        x=rand();
    end
    ch(1) = x;
    for ant=1:(K)
        ch(ant+1)=4*ch(ant)*(1-ch(ant));  %ch相当于βi
        CH(ant,:) = lb+ch(ant)*(ub-lb);    %ub大
        V = (1-setCan)* Alpha_pos+setCan*CH(ant);%V是混沌局部过后生成的空间候选解
        %                 V = Alpha_pos*CH(ant);
        %% 边界控制
        Flag4ub=V>ub;
        Flag4lb=V<lb;
        V=(V.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %% 边界控制结束
         FEs=FEs+1;
        % ObjValV=feval(objfun,V);         %计算函数值
        [FitnessV]=fobj(V);%计算适应度值
        %                 [FitnessV,~,~]=fobj(Positions(i,:)',funcNum);
        
        %    if (ObjValV<GlobalMin)
        if (FitnessV<Alpha_score)
            Alpha_score = FitnessV;
            Alpha_pos = V;
            break;
        end
    end
        
    Convergence_curve(it)=Alpha_score;
    it=it+1; 
end




