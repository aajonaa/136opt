% -------------------------------------------------------------------------- %
%                          适者生存机制 + 差分进化算法                         %
%                                   单纯形法                                  %
%Hybrid Grey Wolf Optimizer Using Elite Opposition-BasedLearning Strategy and Simplex Method
% -------------------------------------------------------------------------- %
function [Leader_pos,Convergence_curve]=WOA_SM_DE(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems
Loser_score = -inf;
Loser_pos = zeros(1,dim);
Convergence_curve=[];    %收敛曲线
FEs = 0;
t = 0;

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
AllFitness = zeros(1,SearchAgents_no);
while FEs<Max_iter
    for i=1:size(Positions,1)%Position数组的行数，其实就是SearchAgent_no种群数量
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;  
        Flag4lb=Positions(i,:)<lb;   
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        FEs = FEs + 1;
        AllFitness(i) = fitness;
        % Update the leader
        if fitness<Leader_score 
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end
        if fitness > Loser_score
            Loser_score = fitness;
%           Loser_pos = Positions(i,:);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [sortFitness,index] = sort(AllFitness,'descend');
    Secondary_score = sortFitness(SearchAgents_no - 1);
    Secondary_pos = Positions(index(SearchAgents_no - 1),:);
    tmpPositions = zeros(SearchAgents_no,dim);
    for i = 1:SearchAgents_no
        tmpPositions(i,:) = Positions(index(i),:);
    end
    Positions = tmpPositions;
    re=rand();   
    reNum = floor(SearchAgents_no * re + 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%引入差分算法的交叉变异，选择机制%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = rand() * (1 - FEs / Max_iter)^2;
    for i=1:reNum
        tmpP(i,:) = Positions(i,:);
        r = randperm(SearchAgents_no,3);
        V = Positions(r(1),:) + rand() * (Positions(r(2),:) - Positions(r(3),:));
        [DEPositions(i,:),Leader_pos,Leader_score,Secondary_pos,Secondary_score,Loser_score,FEs] = DEfun(Positions(i,:),V,Leader_pos,Leader_score,Secondary_pos,Secondary_score,Loser_score,f,dim,FEs,fobj);
    end
    DEP = [DEPositions;tmpP];
    for i =1:size(DEP,1)
        DEPfitness(i) = fobj(DEP(i,:));
        FEs = FEs + 1;        
    end
    [~,index] = sort(DEPfitness);
    Positions(1:reNum,:) = DEP(index(1:reNum),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    a = 2-FEs*((2)/Max_iter);
    a2=-1+FEs*((-1)/Max_iter);

    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        %计算A与C的时候为什么用的不是同一个随机数呢？
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        %w = 1/5 * (3*t^3/500^3 + 2*t^2/500^2);
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)，l是螺旋方式更新位置式中的一个参数，介于[-1,1]之间

        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2) %矩阵position的列数  
            
            if p<0.5   
                if abs(A)>=1       
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    %D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);     
            end        
        end
       %%%%%%%%%%%%%%%%%%%利用单纯形法更新位置向量%%%%%%%%%%%%%%%%%%%%%%
        Pfitness = fobj(Positions(i,:));
        FEs = FEs + 1;
        [SimplexPosition,Leader_score,Leader_pos,FEs] = Simplex(Leader_pos,Leader_score,Secondary_pos,Pfitness,Positions(i,:),FEs,fobj);
        SimplexFitness = fobj(SimplexPosition);
        FEs = FEs + 1;
        if SimplexFitness < Leader_score
            Leader_score = SimplexFitness;
            Leader_pos = SimplexPosition;
        end
        if SimplexFitness > Loser_score
            Loser_score = SimplexFitness;
        end 
        if Pfitness < Leader_score
            Leader_pos = Positions(i,:);
            Leader_score = Pfitness;
        end
        if Pfitness > Loser_score
            Loser_score = Pfitness;
        end
        if SimplexFitness < Pfitness
            Positions(i,:) = SimplexPosition;
        end
        if SimplexFitness > Loser_score
            Loser_score = SimplexFitness;
            %Loser_pos = SimplexPosition;
        end
        if SimplexFitness > Leader_score && SimplexFitness < Secondary_score
            Secondary_pos = SimplexPosition;
            Secondary_score = SimplexFitness;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    t = t + 1;
    Convergence_curve(t) = Leader_score;
end
end
function [Loser_pos,Leader_score,Leader_pos,FEs] = Simplex(Leader_pos,Leader_score,Secondary_pos,Loser_score,Loser_pos,FEs,fobj)
    alpha = 1;
    gamma = 2;
    beta = 0.5;
    Xc = (Leader_pos + Secondary_pos) / 2;
    Xr = Xc + alpha * (Xc - Loser_pos);
    XRfitness = fobj(Xr);
    FEs = FEs + 1;
    if XRfitness < Leader_score
        Leader_score = XRfitness;
        Leader_pos = Xr;
        Xe = Xc + gamma * (Xr - Xc);
        XeFitness = fobj(Xe);
        FEs = FEs + 1;
        if XeFitness < Leader_score
            Leader_score = XeFitness;
            Leader_pos = Xe;
            Loser_pos = Xe;
        else
            Loser_pos = Xr;
        end
    end
    if XRfitness > Loser_score
        Xt = Xc + beta * (Loser_pos - Xc);
        XtFitness = fobj(Xt);
        FEs = FEs + 1;   
        if XtFitness < Leader_score
            Leader_score = XtFitness;
            Leader_pos = Xt;
            Loser_pos = Xt;
        else
            Loser_pos = Xr;
        end
    end
    if Leader_score < XRfitness && XRfitness < Loser_score
        Xw = Xc - beta * (Loser_pos - Xc);
        XwFitness = fobj(Xw);
        FEs = FEs + 1;
        if XwFitness < Leader_score
            Leader_score = XwFitness;
            Leader_pos = Xw;
        end
        if XwFitness < Loser_score
            Loser_pos = Xw;
        else
            Loser_pos = Xr;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%差分算法%%%%%%%%%%%%%%%%%%%%%%%%%
function [VPosition,Leader_pos,Leader_score,Secondary_pos,Secondary_score,Loser_score,FEs] = DEfun(Position,V,Leader_pos,Leader_score,Secondary_pos,Secondary_score,Loser_score,f,dim,FEs,fobj)
    CR = 0.2;
    j = randperm(dim,1);
    for ant = 1:dim
        if rand() <= CR || ant == j
            Position(ant) = V(ant);
        else
            Position(ant) = Position(ant);
        end
    end
    r = rand();
    if r >= 0.5
        VPosition = Position + f * (Position - Leader_pos);
    else
        VPosition = Position + f * (Leader_pos - Position);
    end
    Pfitness = fobj(Position);
    FEs = FEs + 1;
    if Pfitness < Leader_score
        Leader_pos = Position;
        Leader_score = Pfitness;
    end
    if Pfitness > Leader_score && Pfitness < Secondary_score
        Secondary_score = Pfitness;
        Secondary_pos = Position;
    end
    if Pfitness > Loser_score
        Loser_score = Pfitness;
    end
    Vfitness = fobj(VPosition);
    FEs = FEs + 1;
    if Vfitness > Leader_score && Vfitness < Secondary_score
        Secondary_score = Vfitness;
        Secondary_pos = VPosition;
    end
    if Vfitness > Loser_score
        Loser_score = VPosition;
    end
    if Vfitness < Leader_score
        Leader_score = Vfitness;
        Leader_pos = VPosition;
    end
    if Vfitness < Position
        VPosition = VPosition;
    else
        VPosition = Position;
    end
end





