%53框架下 初始化，适应度评价，排序，饥饿更新，权重更新，位置更新
function [bestPositions,Convergence_curve]=HGS(N,Max_iter,lb,ub,dim,fobj)
% disp('HGS 53 is now tackling your problem')

tic
% initialize position
bestPositions=zeros(1,dim);
tempPosition=zeros(N,dim);

Destination_fitness=inf;%change this to -inf for maximization problems
Worstest_fitness=-inf;
AllFitness = inf*ones(N,1);%record the fitness of all positions
VC1 = ones(N,1);%record the variation control of all positions

weight3 = ones(N,dim);%hungry weight of each position
weight4 = ones(N,dim);%hungry weight of each position

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
Convergence_curve=[];

it=1; %Number of iterations
FEs=0;   %Number of current evaluations

hungry = zeros(1,size(X,1));%record the hungry of all positions
count=0;

% Main loop
while  FEs < Max_iter % 3000
    VC2 = 0.03; %The variable of variation control 

    sumHungry = 0;%record the sum of each hungry 
    
    %sort the fitness
    for i=1:size(X,1) % 30
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));  % fobj0
        FEs = FEs+1;

    end
    
    
    
    
    
    [AllFitnessSorted,IndexSorted] = sort(AllFitness);
    bestFitness = AllFitnessSorted(1);
    worstFitness = AllFitnessSorted(size(X,1));
   
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(IndexSorted(1),:);
        Destination_fitness = bestFitness;
        count=0;
    end
    
    if worstFitness > Worstest_fitness
        Worstest_fitness = worstFitness;
    end
    
    
    
    
    
    for i = 1:size(X,1)   %更新饥饿值
         %calculate the variation control of all positions
         VC1(i) = sech(abs(AllFitness(i)-Destination_fitness));    
         %calculate the hungry of each position
        if Destination_fitness == AllFitness(i)          %公式2.7
            hungry(1,i) = 0;
            count = count+1;
            tempPosition(count,:)=X(i,:);
        else                                              %公式2.7
            temprand = rand();
            c = (AllFitness(i)-Destination_fitness)/(Worstest_fitness-Destination_fitness)*temprand*2*(ub-lb);
            if c<100
                b=100*(1+temprand);
            else
                b=c;
            end   
            hungry(1,i) = hungry(1,i)+ max(b); 
            sumHungry = sumHungry + hungry(1,i);
        end
    end 
    
    
    
    
    
    
    %calculate the hungry weight of each position            更新权重
    for i=1:size(X,1)%在这一部分中，用所提出的模型从数学上模拟了个体在搜索过程中的饥饿特征。
        for j=2:size(X,2)
                weight3(i,j) = (1-exp(-abs(hungry(1,i)-sumHungry)))*rand()*2;%公式2.6
                if rand()<VC2
                    weight4(i,j) = hungry(1,i)*size(X,1)/sumHungry*rand();%公式2.5
                else
                    weight4(i,j) = 1;
                end
        end
        
    end
    
    
    
    
    
    % Update the Position of search agents         更新每个个体的位置
    shrink=2*(1-FEs/Max_iter); % a decreases linearly fron 2 to 0
    for i=1:size(X,1)
        if rand<VC2                %公式2.1位置更新策略
            X(i,:) = X(i,j)*(1+randn(1));%表示agent如何在当前位置饥饿且随机地搜索食物;     2.1
        else
            A = randi([1,count]);
            for j=1:size(X,2)
                r = rand();
                vb = 2*shrink*r-shrink;%[-a,a]
                % Moving based on the bestPosition
                % The transformation range is controlled by weight3,bestPositions and X
                if r>VC1(i)
                    X(i,j) = weight4(i,j)*tempPosition(A,j)+vb*weight3(i,j)*abs(tempPosition(A,j)-X(i,j));%2.1
                else
                    X(i,j) = weight4(i,j)*tempPosition(A,j)-vb*weight3(i,j)*abs(tempPosition(A,j)-X(i,j));%2.1
                end
            end
        end
    end

    Convergence_curve(it)=Destination_fitness;
    it=it+1;  % 100 — 50
end
toc
end



