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
function [Alpha_pos,Convergence_curve]=UFGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)% 返回值已经修改成两个

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

CompressOption = 0;

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
% 
% Convergence_curve=zeros(1,Max_iter);%收敛曲线

Convergence_curve=[];%收敛曲线
l=0;% Loop counter
Fes = 0;


NumP = size(Positions,1);
Father = 1:NumP;
MST1 = [];
MST2 = [];
MST3 = [];

Pos_Temp = initialization(SearchAgents_no,dim,ub,lb);
% for i=1:size(Positions,1)
%     for j=1:size(Positions,2)
%         r = exp(j);
%         Pos_Temp(i,j) = i * r - floor(i * r);
%         Positions(i,j) =lb + Pos_Temp(i,j) * (ub-lb);%
%     end
% end
  
FitnessVal = zeros(1,SearchAgents_no);
    

% Main loop
while Fes<Max_iter
   
    for i=1:size(Positions,1)
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        if Fes < Max_iter
            Fes = Fes + 1;%每次评估要使评估值加一
            % Calculate objective function for each search agent
            fitness = fobj(Positions(i,:));
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
%    for i=1:SearchAgents_no
%        FitnessVal = fobj(Positions);
%     end
%     Destination_fitness=get_fitness(X,fobj);
%     [Destination_fitness_sorted,n]=sort(Destination_fitness);
%     sorted_X=X(n,:);     
    X = Positions;
   
    [FitnessVal,Fes] = get_fitness(X,fobj,Fes);
    [~,index] = sort(FitnessVal);
    Positions = X(1:index,:);
   
    %%
    MST1 = [MST1;Alpha_pos];
    MST2 = [MST2;Beta_pos];
    MST3 = [MST3;Delta_pos];
    for i=1:size(Positions,1)
        Fes = Fes + 1;
        fitness_temp = fobj(Positions(i,:));
        Pos_temp = Positions(i,:);
        
        if fitness_temp > Delta_score
            [Node1,Node2]=GetSquarePos(index(i),NumP);
            if mod(i,3) == 0
                Father = MergeNode(Node1,Node2,Father,CompressOption);
                MST1 = [MST1;Pos_temp];
            end
            if mod(i,3) == 1
                %[Node1,Node2]= GetSquarePos(index(i),NumP);
                Father =  MergeNode(Node1,Node2,Father,CompressOption);
                MST2 = [MST2;Pos_temp ];
            end
            if mod(i,3) == 2
                %[Node1,Node2]= GetSquarePos(index(i),NumP);
                Father = MergeNode(Node1,Node2,Father,CompressOption);
                MST3 = [MST3;Pos_temp ];
            end
        end
        
    end
    %
    if size(MST1,1)>SearchAgents_no
        [FitnessM1,Fes] = get_fitness(MST1,fobj,Fes);
        [~,indexM1] = sort(FitnessM1);
        MST1 = MST1(indexM1,:);
        MST1 = MST1(1:SearchAgents_no,:);%将 Positions按从小到大的顺序排序
    end
    
    if size(MST2,1)>SearchAgents_no
        [FitnessM2,Fes] = get_fitness(MST2,fobj,Fes);
        [~,indexM2] = sort(FitnessM2);
        MST2 = MST2(indexM2,:);
        MST2 = MST2(1:SearchAgents_no,:);%将 Positions按从小到大的顺序排序
    end
    
    if size(MST3,1)>SearchAgents_no
       [FitnessM3,Fes] = get_fitness(MST3,fobj,Fes);
        [~,indexM3] = sort(FitnessM3);
        MST3 = MST3(indexM3,:);
        MST3 = MST3(1:SearchAgents_no,:);%将 Positions按从小到大的顺序排序
    end
    %
 
    a=2-Fes*((2)/Max_iter); % a decreases linearly fron 2 to 0
    %a = 2 * sqrt(1-(l/2)^2);
    % Update the Position of search agents including omegas
    [MST1,Fes] = update_w(Alpha_pos,Beta_pos,Delta_pos,a,fobj,MST1,Fes);
	[MST2,Fes]  = update_w(Alpha_pos,Beta_pos,Delta_pos,a,fobj,MST2,Fes);
	[MST3,Fes] = update_w(Alpha_pos,Beta_pos,Delta_pos,a,fobj,MST3,Fes);
    fitnessF = zeros(1,size(MST1,1));
    fitnessS = zeros(1,size(MST2,1));
    fitnessH = zeros(1,size(MST3,1));
    for i=1:size(MST1,1)
        Fes=Fes+1;
        fitnessF(i) = fobj(MST1(i,:));
    end
    for i=1:size(MST2,1)
        Fes=Fes+1;
        fitnessS(i) = fobj(MST2(i,:));
    end
    for i=1:size(MST3,1)
        Fes=Fes+1;
        fitnessH(i) = fobj(MST3(i,:));
    end







%     [~,index1] = sort(fitnessF);
%     [~,index2] = sort(fitnessS);
%     [~,index3] = sort(fitnessH);
    Positions_temp = [Positions;MST1;MST2;MST3];
      % for i=1:SearchAgents_no
      
     [FitnessVal,Fes] = get_fitness(Positions_temp,fobj,Fes);
    [~,index] = sort(FitnessVal);
    Positions_temp = Positions_temp(index,:);
    if index > SearchAgents_no
        Positions = Positions_temp(1:SearchAgents_no,:);
    else
        Positions = Positions_temp(1:index,:);
    end
    
    
%     FitnessVal = fobj(Positions_temp);
%     [~,indexx] = sort(FitnessVal);
%     Positions = Positions_temp(indexx,:);
%    Positions = Positions(1:SearchAgents_no,:);%将 Positions按从小到大的顺序排序
%    f1 = min(fitnessF);
%    f2 = min(fitnessS);
%    f3 = min(fitnessH);
%    Alpha_score = min(f1,min(f2,f3));
    
    l=l+1;
    Convergence_curve(l) = Alpha_score;
end
end


function [Y,Fes]=get_fitness(X,fobj,Fes)
y=zeros(1,size(X,1));
for i=1:size(X,1)
      Fes = Fes + 1;%每次评估要使评估值加一
    y(1,i)=fobj(X(i,:));
end
Y=y;
end


function Node = RootNode(Node,Father)  % 找到根结点
    %Father = zeros(SearchAgents_no,dim);
    temp = size(Father);
   %x = temp(2)
    for i=1:size(Node)
       % Node(i)
        if Node(i) > temp(2)
            Node(i) = mod(Node(i),temp(2));
            if Node(i) == 0
                Node(i) = 1;
            end
           % Node(i) =  randperm(temp(2),1);
        end
        %Node(i)
        while(Node(i) ~= Father(Node(i)))
            Node(i) = Father(Node(i));
        end
    end
    
end


function Father = MergeNode(Node1,Node2,Father,CompressOption)
   % Father =[];
    if nargin == 3
        CompressOption = 1;
    end
    RootNode1 = RootNode(Node1,Father);
    if CompressOption  % 采取压缩操作
        while(Node1~=RootNode1)  % 压缩Node1所在的集合
            t = Father(Node1);
            Father(Node1) = RootNode1;
            Node1 = t;
        end
        while(Node2~=Father(Node2))    % 压缩Node2所在的集合
            t = Father(Node2);
            Father(Node2) = RootNode1; % 全部改变为Node1的根结点
            Node2 = t;
        end
        Father(Node2) = RootNode1; % 改变根结点
    else
        Father(RootNode(Node2,Father)) = RootNode1; % 改变根结点
    end
end


function [Node1,Node2] = GetSquarePos(pos,row )  % 输入向量中的位置，输出矩阵中的位置
    %S = cumsum([0:Dim-1])-pos;
    %Node1 = min(find(S>=0));         % 列数
    %Node2 = pos-sum(0:Node1-2);      % 行数
    Node1 = floor(pos/row);
    if Node1 == 0
        Node1 = 1;
    end
    Node2 = mod(pos,row);
    if Node2 == 0
        Node2 = 1;
    end
end


function [X,Fes] = update_w(Alpha_pos,Beta_pos,Delta_pos,a,fobj,MST,Fes)
for i=1:size(MST,1)
    for j=1:size(MST,2)
        
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A1=2*a*r1-a; % Equation (3.3)
        C1=2*r2; % Equation (3.4)
        
        D_alpha=abs(C1*Alpha_pos(j)-MST(i,j)); % Equation (3.5)-part 1
        X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
        
        r1=rand();
        r2=rand();
        
        A2=2*a*r1-a; % Equation (3.3)
        C2=2*r2; % Equation (3.4)
        
        D_beta=abs(C2*Beta_pos(j)-MST(i,j)); % Equation (3.5)-part 2
        X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2
        
        r1=rand();
        r2=rand();
        
        A3=2*a*r1-a; % Equation (3.3)
        C3=2*r2; % Equation (3.4)
        
        D_delta=abs(C3*Delta_pos(j)-MST(i,j)); % Equation (3.5)-part 3
        X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3
        
        %Positions(i,j)=(X1+X2*0.67+X3*0.33)/3;% Equation (3.7)
        
        w1 = abs(X1)/abs(X1+X2+X3);
        w2 = abs(X2)/abs(X1+X2+X3);
        w3 = abs(X3)/abs(X1+X2+X3);
        MST(i,j)=(X1*w1+X2*w2+X3*w3)/3;% Equation (3.7)
        
         Fes=Fes+2;%每次评估要使评估值加一
        Moth_pos_m_gaus=MST(i,:)*(1+randn(1));%经过高斯变异后搜索代理的位置
        Moth_fitness_m_gaus=fobj(Moth_pos_m_gaus);%经过高斯变异后搜索代理的适应度值
        Moth_fitness_s=fobj(MST(i,:));%没有经过高斯变异的搜索代理的适应度值
        
        Moth_fitness_comb=[Moth_fitness_m_gaus,Moth_fitness_s];%将两次的适应度值拼接在一起
        [~,m]=min(Moth_fitness_comb);%找到所有适应度中的最小值，m为最小值的下标
        if m==1%如果m为1则表示经过高斯变异后的位置比较好，用高斯变异后的位置取代原来的位置，否则就原来的位置不变
            MST(i,:)=Moth_pos_m_gaus;
        else
            MST(i,:)=MST(i,:);
        end   
     
    end
        
end
X = MST;
end


