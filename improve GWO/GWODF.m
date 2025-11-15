

% Grey Wolf Optimizer
function [Alpha_pos,Convergence_curve]=GWODF(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

it=1;FEs=0;
Convergence_curve=[];
Maximum_Diffusion = 2;
Walk = 0.75;
% Main loop
while  FEs < MaxFEs
% % Main loop
% while l<Max_iter
    for i=1:size(Positions,1)  
        
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        if FEs<MaxFEs
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
    
    
    a=2-FEs*((2)/MaxFEs); % a decreases linearly fron 2 to 0
    New_Positions = [];%creating new Positions
    FitVector = [];%creating vector of fitness functions
    %diffusion process occurs for all Positionss in the group
    for i = 1 : size(Positions,1)
        %creating new Positionss based on diffusion process
        [NP, fit] = Diffusion_Process(Positions(i,:),Maximum_Diffusion,Walk,lb,ub,FEs,Alpha_pos,fobj);
        New_Positions = [New_Positions;NP];
        FitVector = [FitVector,fit];
    end
    Positions = New_Positions;
    % Update the Position of search agents including omegas
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
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
    end
    
     %%  %使用搜索策略预形成混沌局部搜索
    setCan = (MaxFEs-FEs+1)/MaxFEs;%相当于入
    x = rand();
    while(~(x~=0.25 && x~=0.5 && x~=0.75 && x~=1))
        x=rand();
    end
    ch(1) = x;
    K = SearchAgents_no;
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
end



function [createPosition, fitness] = Diffusion_Process(Position,Maximum_Diffusion,Walk,lb,ub,g,BestPosition,fobj)
    %calculating the maximum diffusion for each Position
    NumDiffiusion = Maximum_Diffusion;
    New_Position = Position;
    
    %Diffiusing Part*******************************************************
    for i = 1 : NumDiffiusion
        %consider which walks should be selected.
        if rand < Walk 
            GeneratePosition = normrnd(BestPosition, (log(g)/g)*(abs((Position - BestPosition))), [1 size(Position,2)]) + ...
                (randn*BestPosition - randn*Position);
        else
            GeneratePosition = normrnd(Position, (log(g)/g)*(abs((Position - BestPosition))),...
                [1 size(Position,2)]);
        end
        New_Position = [New_Position;GeneratePosition];
    end
    %check bounds of New Position
    New_Position = Bound_Checking(New_Position,lb,ub);
    %sorting fitness
    fitness = [];
    for i = 1 : size(New_Position,1)
        g = g + 1;
        fitness = [fitness;fobj(New_Position(i,:))];
    end
    [fit_value,fit_index] = sort(fitness);
    fitness = fit_value(1,1);
    New_Position = New_Position(fit_index,:);
    createPosition = New_Position(1,:);
    %======================================================================
end

function V = Bound_Checking(p,lb,ub)
%     for i = 1 : size(p,1)
%         upper = double(gt(p(i,:),upB));
%         lower = double(lt(p(i,:),lowB));
%         up = find(upper == 1);
%         lo = find(lower == 1);
%         if (size(up,2)+ size(lo,2) > 0 )
%             for j = 1 : size(up,2)
%                 p(i, up(j)) = (upB(up(j)) - lowB(up(j)))*rand()...
%                     + lowB(up(j));
%             end
%             for j = 1 : size(lo,2)
%                 p(i, lo(j)) = (upB(lo(j)) - lowB(lo(j)))*rand()...
%                     + lowB(lo(j));
%             end
%         end
%     end
for i=1:size(p,1)
    Flag4ub=p(i,:)>ub;
    Flag4lb=p(i,:)<lb;
    p(i,:)=(p(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
end
V = p;
end