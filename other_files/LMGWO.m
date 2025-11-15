
% Grey Wolf Optimizer
function [Alpha_pos,Convergence_curve]=LMGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=[];

% Loop counter
ML = 10;
MemoryLimit = ML;
FEs = 0;
iter = 0;
gBestCollPos=zeros(ML,dim);
gBestCollCost=1.0e+100*ones(ML,1);
gBestCollProb=zeros(ML,1);
gBestCollCount=0;

while FEs<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
       % display(gBestCollCont);
        % Pos = Positions;
        % Calculate objective function for each search agent
        if FEs < Max_iter
            FEs = FEs + 1;
            fitness=fobj(Positions(i,:));
            fitness_val(i) = fitness;
            %fitnessVal(i) = fobj(Positions(i,:));
            %[FitnessVal_oder,index] = sort(fitnessVal);%FitnessHe2保存的是排序后的适应度值，index是FitnessHe2中保存的FitnessHe1的索引
            % Positions = Pos(index,:);
            % Positions = Positions(1:SearchAgents_no,:);%将 Positions按从小到大的顺序排序

            % Update Alpha, Beta, and Delta
            if fitness<Alpha_score
                Alpha_score=fitness; % Update alpha
                Alpha_pos=Positions(i,:);
               % display(gBestCollCont);
                if gBestCollCount <= ML && gBestCollCount>0
                    gBestCollPos(gBestCollCount,:) = Positions(i,:);
                    gBestCollCost(gBestCollCount) = fitness;
                     %display(gBestCollCont);
                else
                   % display(gBestCollCont);
                    gBestCollCount = 1;
                    gBestCollPos(gBestCollCount,:) = Positions(i,:);
                    gBestCollCost(gBestCollCount) = fitness;
                end
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
    F=zeros(ML,1);
    MeanCost = mean(gBestCollCost);
    for j=1:MemoryLimit
        F(j) = exp(-gBestCollCost(j)/MeanCost);
    end
    gBestCollProb(:)=F/sum(F);
    
    if range(gBestCollCost)~=0
        pr=RouletteWheelSelection(gBestCollProb);
        % Choose k randomly, not equal to i
        K=[1:pr-1 pr+1:MemoryLimit];
        k=K(randi([1 numel(K)]));
        gbestX=gBestCollPos(k,:);
    else
        gbestX = Alpha_pos;
    end
    
    %  a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    n=1-FEs*((1)/Max_iter);
    if FEs <= Max_iter/2
        a = 0 + (2 - 0)*(1 + (cos(((FEs-1)*pi)/(Max_iter-1)))^n)/2;
    else
        a = 0 + (2 - 0)*(1 - abs(cos(((FEs-1)*pi)/(Max_iter-1)))^n)/2;
    end
    
    %% 
    
    
    %%
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)
            
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*gbestX(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=gbestX(j)-A1*D_alpha; % Equation (3.6)-part 1
            
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
            
            %Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            w1 = abs(X1)/abs(X1+X2+X3);
            w2 = abs(X2)/abs(X1+X2+X3);
            w3 = abs(X3)/abs(X1+X2+X3);
            Positions(i,j)=(X1*w1+X2*w2+X3*w3)/3;% Equation (3.7)
        end
    end
    iter = iter + 1;
    Convergence_curve(iter)=Alpha_score;
end
end


