% The Whale Optimization Algorithm
function [Leader_pos,Convergence_curve]=RDWOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)


K = SearchAgents_no;
% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
% Convergence_curve=zeros(1,MaxFEs);
%Convergence_curve(0)=0;
Convergence_curve=[];
FEs=0;
t=0;
s=0;

% Main loop
while  FEs < MaxFEs
 
    for i=1:size(Positions,1)    

        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));
        FEs=FEs+1;
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
            s=s/2;           % s作为陷入局部最优程度的一个间接表示量
        end
        s=s+1;

    end

    
    %% 随机替换，在算法后期将当前位置以一定概率替换成最优解的位置，然后进行评估
    for l=1:size(Positions,1)
        M=Positions(l,:);
%        if(FEs/MaxFEs>0.7 && FEs <0.9*MaxFEs)      %在算法进行到70%-90%阶段有概率替换
          for h=1:size(Positions,2)
              if(tan(pi*(rand-0.5))<(1-FEs/MaxFEs))  %根据算法剩余运行次数占总运行次数的比值与柯西随机数相比较，使当前位置有一定几率向最优位置靠拢，越后期替换概率越小
                    M(h)= Leader_pos(h);               %%    
              end
          end   
                Fitnessm=fobj(M);               %计算适应度值
                  FEs=FEs+1;
            if (Fitnessm<Leader_score)
               Leader_score = Fitnessm;
               Leader_pos =M;
               break;
            end      
%             Convergence_curve(FEs)=Leader_score;
%        else
%            continue;
%        end
    end
 

   %% 加入双自适应权重，w使算法在前期有较好的全局优化能力，w1使后期有较好的搜索能力
    w1=(1-FEs/MaxFEs)^(1-tan(pi*(rand-0.5))*s/MaxFEs);   %权重w会根据算法陷入局部最优程度呈曲线递减
    w2=(2-2*FEs/MaxFEs)^(1-tan(pi*(rand-0.5))*s/MaxFEs);  %权重w1会根据算法陷入局部最优程度呈曲线递减
     a=2-FEs*((2)/MaxFEs); % a decreases linearly fron 2 to 0 in Eq. (2.3)

    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+FEs*((-1)/MaxFEs);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    if(FEs/MaxFEs>0.5)                                                          %在算法后半段进行更小范围的局部搜索
                    Positions(i,j)=X_rand(j)-A*D_X_rand*w2;      % Eq. (2.8)
                    else
                    Positions(i,j)=X_rand(j)*w1-A*D_X_rand;      % Eq. (2.8)                      %在算法前半段进行更大范围的全局搜索
                    end
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                        if(FEs/MaxFEs>0.5)                                                       %在算法后半段进行更小范围的局部搜索
                    Positions(i,j)=Leader_pos(j)-w2*A*D_Leader;      % Eq. (2.2)
                        else
                          Positions(i,j)=w1*Leader_pos(j)-A*D_Leader;      % Eq. (2.2)            %在算法前半段进行更大范围的全局搜索
                        end  
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                if(FEs/MaxFEs>0.5)                                                               %在算法后半段进行更小范围的局部搜索
                Positions(i,j)=w2*distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                else
                   Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+w1*Leader_pos(j);       %在算法前半段进行更大范围的全局搜索
                end
                
            end
            
        end
    end
    t=t+1;
     Convergence_curve(t)=Leader_score;
    

end



