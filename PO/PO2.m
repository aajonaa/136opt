%53框架PO
%_________________________________________________________________________________
%  Political Optimizer: A novel socio-inspired meta-heuristic 
%                       for global optimization source codes version 1.0
%
%  Developed in MATLAB R2015a
%
%  Author and programmer: Qamar Askari
%
%         e-Mail: l165502@lhr.nu.edu.pk
%                 syedqamar@gift.edu.pk
%
%
%   Main paper:
%   Askari, Q., Younas, I., & Saeed, M. (2020). Political Optimizer: 
%       A novel socio-inspired meta-heuristic for global optimization.
%   Knowledge-Based Systems, 2020, 
%   DOI: https://doi.org/10.1016/j.knosys.2020.105709
%____________________________________________________________________________________

function [Leader_pos,Convergence_curve]=PO(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
% initialize position vector and score for the leader初始化领导的位置向量和分数
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems
areas=5;
parties=6;         %Attention!!!!!!  dim=areas*parties
lambda=1;
FES=0;
%Initialize the positions of search agents初始化搜索代理的位置
Positions=initialization(SearchAgents_no,dim,ub,lb);
auxPositions = Positions;
prevPositions = Positions;
Convergence_curve=[];
fitness=zeros(SearchAgents_no, 1);

%Running phases for initializations初始化的运行阶段
%Election选举   %Run election phase选举运行阶段
for i=1:size(Positions,1)
	% Return back the search agents that go beyond the boundaries of the search space返回超出搜索空间边界的搜索代理
 	Flag4ub=Positions(i,:)>ub;
 	Flag4lb=Positions(i,:)<lb;
 	Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
       
	%Calculate objective function for each search agent计算每个搜索代理的目标函数
	fitness(i,1)=fobj(Positions(i,:));
    FES=FES+1;    
	%Update the leader更新领导者
    if fitness(i,1)<Leader_score % Change this to > for maximization problem
        Leader_score=fitness(i,1); 
        Leader_pos=Positions(i,:);
    end        
end
auxFitness = fitness;
prevFitness = fitness;
%Government Formation                        一。政党形成
aWinnerInd=zeros(areas,1);   %Indices of area winners in x  x区域优胜者指数
aWinners = zeros(areas,dim); %Area winners are stored separately区域优胜者分开存放
for a = 1:areas
	[aWinnerFitness,aWinnerParty]=min(fitness(a:areas:SearchAgents_no));
	aWinnerInd(a,1) = (aWinnerParty-1) * areas + a;
    aWinners(a,:) = Positions(aWinnerInd(a,1),:);
end    

%Finding party leaders
pLeaderInd=zeros(parties,1);    %Indices of party leaders in x
pLeaders = zeros(parties,dim);  %Positions of party leaders in x
for p = 1:parties
	pStIndex = (p-1) * areas + 1;
	pEndIndex = pStIndex + areas - 1;
	[partyLeaderFitness,leadIndex]=min(fitness(pStIndex:pEndIndex)); 
	pLeaderInd(p,1) = (pStIndex - 1) + leadIndex; %Index of party leader
    pLeaders(p,:) = Positions(pLeaderInd(p,1),:);
end
t=0;% Loop counter循环计数器


while FES<MaxFEs
    prevFitness = auxFitness;
    prevPositions = auxPositions;
    auxFitness = fitness;
    auxPositions = Positions;  %对于政党成员们的当前位置和适应度的暂时复制
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

    %Election Campaign                         二。竞选活动    在这个阶段中更新政党成员的位置
    for whichMethod = 1:2
        for a = 1:areas
            for p = 1:parties
                i = (p-1)*areas + a; %index of member成员索引
        
                for j=1:dim
                    if whichMethod == 1         %position-updating w.r.t party leader职位更新w.r.t党魁
                        center = pLeaders(p,j);
                    elseif whichMethod == 2     %position-updating w.r.t area winner位置更新w.r.t区域赢家
                        center = aWinners(a,j);
                    end
            
                %Cases of Eq. 9 in paper
                    if prevFitness(i) >= fitness(i) 
                        if (prevPositions(i,j) <= Positions(i,j) && Positions(i,j) <= center) ...
                                || (prevPositions(i,j) >= Positions(i,j) && Positions(i,j) >= center)
                    
                            radius = center - Positions(i,j); 
                            Positions(i,j) = center + rand() * radius;
                        elseif (prevPositions(i,j) <= Positions(i,j) && Positions(i,j) >= center && center >= prevPositions(i,j)) ...
                                || (prevPositions(i,j) >= Positions(i,j) && Positions(i,j) <= center && center <= prevPositions(i,j))
                    
                            radius = abs(Positions(i,j) - center);
                            Positions(i,j) = center + (2*rand()-1) * radius;
                        elseif (prevPositions(i,j) <= Positions(i,j) && Positions(i,j) >= center && center <= prevPositions(i,j)) ...
                                || (prevPositions(i,j) >= Positions(i,j) && Positions(i,j) <= center && center >= prevPositions(i,j))
                    
                            radius = abs(prevPositions(i,j) - center);
                            Positions(i,j) = center + (2*rand()-1) * radius;
                        end
                
                %Cases of Eq. 10 in paper
                   elseif prevFitness(i) < fitness(i) 
                        if (prevPositions(i,j) <= Positions(i,j) && Positions(i,j) <= center) ...
                                || (prevPositions(i,j) >= Positions(i,j) && Positions(i,j) >= center)
                    
                    
                            radius = abs(Positions(i,j) - center);
                            Positions(i,j) = center + (2*rand()-1) * radius;
                        elseif (prevPositions(i,j) <= Positions(i,j) && Positions(i,j) >= center && center >= prevPositions(i,j)) ...
                                || (prevPositions(i,j) >= Positions(i,j) && Positions(i,j) <= center && center <= prevPositions(i,j))
                
                            radius = Positions(i,j) - prevPositions(i,j);
                            Positions(i,j) = prevPositions(i,j) + rand() * radius;
                        elseif (prevPositions(i,j) <= Positions(i,j) && Positions(i,j) >= center && center <= prevPositions(i,j)) ...
                                || (prevPositions(i,j) >= Positions(i,j) && Positions(i,j) <= center && center >= prevPositions(i,j))
            
                            center2 = prevPositions(i,j);
                            radius = abs(center - center2);
                            Positions(i,j) = center + (2*rand()-1) * radius;
                        end
                    end
            
                end
            end
        end
    end
    %Party Switching                              三。政党转换
    psr = (1-FES*((1)/MaxFEs)) * lambda;

    for p=1:parties
        for a=1:areas
            fromPInd = (p-1)*areas + a;
            if rand() < psr
                %Selecting a party other than current where want to send the member选择要将成员发送到当前位置以外的参与方
                toParty = randi(parties);
                while(toParty == p)
                    toParty = randi(parties);
                end
            
                %Deciding member in TO party党的决定成员
                toPStInd = (toParty-1) * areas + 1;
                toPEndIndex = toPStInd + areas - 1;
                [~,toPLeastFit] = max(fitness(toPStInd:toPEndIndex));
                toPInd = toPStInd + toPLeastFit-1;
            
            
                %Deciding what to do with member in FROM party and switching决定如何处理从党的成员和转换
                fromPInd = (p-1)*areas + a;
                temp = Positions(toPInd,:);
                Positions(toPInd,:) = Positions(fromPInd);
                Positions(fromPInd,:)=temp;
                
                temp = fitness(toPInd);
                fitness(toPInd) = fitness(fromPInd);
                fitness(fromPInd) = temp;
            end
        end
    end
    %Election                                       四。政党间选举
    for i=1:size(Positions,1)
        % Return back the search agents that go beyond the boundaries of the search space返回超出搜索空间边界的搜索代理
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
       
        %Calculate objective function for each search agent计算每个搜索代理的目标函数
        fitness(i,1)=fobj(Positions(i,:));
        FES=FES+1;
        %Update the leader
        if fitness(i,1)<Leader_score % Change this to > for maximization problem
            Leader_score=fitness(i,1); 
            Leader_pos=Positions(i,:);
        end        
    end
    %Government  Formation
    aWinnerInd=zeros(areas,1);   %Indices of area winners in x
    aWinners = zeros(areas,dim); %Area winners are stored separately
    for a = 1:areas
        [aWinnerFitness,aWinnerParty]=min(fitness(a:areas:SearchAgents_no));
        aWinnerInd(a,1) = (aWinnerParty-1) * areas + a;
        aWinners(a,:) = Positions(aWinnerInd(a,1),:);
    end    

    %Finding party leaders
    pLeaderInd=zeros(parties,1);    %Indices of party leaders in x
    pLeaders = zeros(parties,dim);  %Positions of party leaders in x
    for p = 1:parties
        pStIndex = (p-1) * areas + 1;
        pEndIndex = pStIndex + areas - 1;
        [partyLeaderFitness,leadIndex]=min(fitness(pStIndex:pEndIndex)); 
        pLeaderInd(p,1) = (pStIndex - 1) + leadIndex; %Indexof party leader
        pLeaders(p,:) = Positions(pLeaderInd(p,1),:);
    end
    %Parliamentarism议会主义              五。 parliamentary affairs   议会事务
    for a=1:areas
        newAWinner = aWinners(a,:);
        i = aWinnerInd(a);    
    
        toa = randi(areas);
        while(toa == a)
            toa = randi(areas);
        end
        toAWinner = aWinners(toa,:);
        for j = 1:dim
            distance = abs(toAWinner(1,j) - newAWinner(1,j));
            newAWinner(1,j) = toAWinner(1,j) + (2*rand()-1) * distance;
        end
        newAWFitness=fobj(newAWinner(1,:));
        FES=FES+1;
        %Replace only if improves仅在改进时更换     update previous positions and fitness更新之前的位置和适应度
        if newAWFitness < fitness(i) 
            Positions(i,:) = newAWinner(1,:);
            fitness(i) = newAWFitness;
            aWinners(a,:) = newAWinner(1,:);
        end
    end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    t=t+1;
    Convergence_curve(t)=Leader_score;
  %  [t Leader_score];
end
end

