%% Fata morgana algorithm .Qi A
function [bestPos,cg_curve]=FATA_FEs(N,MaxFEs,lb,ub,dim,fobj)
    % initialize position
     worstInte=0; %Parameters of Eq.(4);Worst population quality
     bestInte=Inf;%Parameters of Eq.(4);Best population quality
    noP=N;%Population size
    arf=0.2;%Eq. (15) reflectance=0.2; reflectance was verified experimentally in section 4.2 of the paper
    best_pos=zeros(1,dim);%Globally optimal individual
    cg_curve=[];%Convergent curve vector
    gBestScore=inf;%change this to -inf for maximization problems
    X=initialization(noP,dim,ub,lb);%Initialize the set of random solutions
    fitness=zeros(noP,1)+inf;
    it=1;%Number of iterations   
    FEs=0;
    lb=ones(1,dim).*lb; % lower boundary 
    ub=ones(1,dim).*ub; % upper boundary
    % Main
    while  FEs < MaxFEs 
        for i=1:size(X,1)     
            Flag4ub=X(i,:)>ub;
            Flag4lb=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            FEs=FEs+1;    
            fitness(i)=fobj(X(i,:));
             %Make greedy selections
            if(gBestScore>fitness(i))
                gBestScore=fitness(i);
                best_pos=X(i,:);
            end
        end
        [Order,Index] = sort(fitness);  
        worstFitness = Order(N); 
        bestFitness = Order(1);
     %% The mirage light filtering principle 
     % Definite integral strategy
     %This formula is based on Figure 5 in the paper to calculate the fitness constant integral mass of the population
     Integral=cumtrapz(Order);
     if Integral(N)>worstInte
         worstInte=Integral(N);
     end
      if Integral(N)<bestInte 
         bestInte =Integral(N);
      end
    IP=(Integral(N)-worstInte)/(bestInte-worstInte+eps);% Eq.(4) population quality factor
     %% Calculation zeroo and Para2
        a = tan(-(FEs/MaxFEs)+1);
        b = 1/tan(-(FEs/MaxFEs)+1);
        %% Population location renewal
         for i=1:size(X,1) 
             zeroo=a*rand(1,dim)-a*rand(1,dim); %Parameters of Eq.(10)
             Para2=b*rand(1,dim)-b*rand(1,dim);%Parameters of Eq.(13)
             p=((fitness(i)-worstFitness))/(gBestScore-worstFitness+eps);% Parameters of Eq.(5) individual quality factor
             %% Eq.(1) 
             if  rand>IP 
                 X(i,:) = (ub-lb).*rand+lb;
             else
                for j=1:dim
                    randIndex=floor(rand*N+1);
                    if rand<p   
                        X(i,j)=best_pos(j)+X(i,j).*zeroo(j);%Light refraction(first phase)  Eq.(8)      
                    else   
                        X(i,j)=X(randIndex,j)+Para2(j).*X(i,j);%Light refraction(second phase)   Eq.(11)    
                        X(i,j)=(0.5*(arf+1).*(lb(j)+ub(j))-arf.*X(i,j));%Light total internal reflection Eq.(14)  
                    end      
                end
            end
         end     
        cg_curve(it)=gBestScore;
        it=it+1;
        bestPos=best_pos;
    end
end