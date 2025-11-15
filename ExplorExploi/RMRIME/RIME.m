
% RIME Algorithm modified by Jona
% 2024-5-8 for balance and diversity.
% function [best_pos,Convergence_curve]=RIME(N,Max_FEs,lb,ub,dim,fobj)
function [average_distance,Div,best_pos,Convergence_curve]=RIME(N,Max_FEs,lb,ub,dim,fobj)
    % disp('RIME is now tackling your problem')
    
    %%
    Div=[];
    average_distance=[];
    iter=1;  %Number of iterations
    %%
    
    % initialize position
    best_pos=zeros(1,dim);
    bestFitness=inf;%change this to -inf for maximization problems
    X=initialization(N,dim,ub,lb);%Initialize the set of random solutions
    
    %% Add the sentence after X is defined
    D= pdist(X);
    average_distance(iter)=(sum(D(1:end)))/(N*(N-1)/2);
    Div(iter)=Divergence_calculation(X,dim,N);
    %%
    
    X2 = X;
    Lb=lb.*ones(1,dim);% lower boundary 
    Ub=ub.*ones(1,dim);% upper boundary
    % it=1;%Number of iterations
    FEs=0;
    iter = 1;
    Convergence_curve=[];
    AllFitness=zeros(1,N);%Initialize the fitness value
    newAllFitness=zeros(1,N);
    W = 5;%Soft-rime parameters, discussed in subsection 4.3.1 of the paper
    %Calculate the fitness value of the initial position
    for i=1:N
        AllFitness(1,i)=fobj(X(i,:));%Calculate the fitness value for each search agent
        FEs=FEs+1;
        %Make greedy selections
        if AllFitness(1,i)<bestFitness
            bestFitness=AllFitness(1,i);
            best_pos=X(i,:);
        end
    end
    % Main loop
    while FEs < Max_FEs

        [sorted_AllFitness, ~] = sort(AllFitness);
        Rolette_index=RouletteWheelSelection(1./sorted_AllFitness);
        if Rolette_index==-1  
            Rolette_index=1;
        end

        F = (rand-0.5)*2*cos((pi*FEs/(Max_FEs/10)))*(1-round(FEs*W/Max_FEs)/W);%Parameters of Eq.(3),(4),(5)
        p = sqrt(FEs/Max_FEs);%Eq.(6)
        newX = X;%Recording new populations
        normalizedAllFitness=normr(AllFitness);%Parameters of Eq.(7)
        for i=1:N
            for j=1:dim
                %Soft-rime search strategy
                r1=rand();
                if r1< p
                    newX(i,j)=best_pos(1,j)+ F *((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
                end
                %Hard-rime puncture mechanism
                r2=rand();
                if r2 < normalizedAllFitness(i)
                    newX(i,j)=best_pos(1,j);%Eq.(7)
                end
            end
        end
        for i=1:N
            %Boundary absorption
            Flag4ub=newX(i,:)>ub;
            Flag4lb=newX(i,:)<lb;
            newX(i,:)=(newX(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            newAllFitness(1,i)=fobj(newX(i,:));
            FEs=FEs+1;
            %Positive greedy selection mechanism
            if newAllFitness(1,i)<AllFitness(1,i)
                AllFitness(1,i) = newAllFitness(1,i);
                X(i,:) = newX(i,:);
                if newAllFitness(1,i)< bestFitness
                   bestFitness=AllFitness(1,i);
                   best_pos=X(i,:);
                end
            end
        end
        
        %% Add the sentence before plot the dot.
        % D= pdist(X,'euclidean');
        D= pdist(X);
        average_distance(iter)=(sum(D(1:end)))/(N*(N-1)/2);
        Div(iter)=Divergence_calculation(X,dim,N);
        %%
        
        Convergence_curve(iter)=bestFitness;
        iter=iter+1;
    end
end

function choice = RouletteWheelSelection(weights)
  accumulation = cumsum(weights);
  p = rand() * accumulation(end);
  chosen_index = -1;
  for index = 1:length(accumulation)
    if (accumulation(index) > p)
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
end
