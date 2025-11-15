function [best_pos,Convergence_curve]=MGOmod(N,MaxFEs,lb,ub,dim,fobj)
%% Initialization
FEs=0;
bestFitness=inf; %change this to -inf for maximization problems
best_pos = zeros(1,dim);
X=initialization(N,dim,ub,lb); %Initialize the set of random solutions
AllFitness = zeros(1,N);
for i=1:N
    AllFitness(i)=fobj(X(i,:)) ;
    FEs=FEs+1;
    if AllFitness(i)<bestFitness 
        best_pos=X(i,:); 
        bestFitness = AllFitness(1,i);
    end
end

Convergence_curve=[];
it=1;
rec = 1;

w = 2;
rec_num = 10;
divide_num = dim/4;
d1 = 0.2;

newX = zeros(N, dim);
newAllFitness = zeros(1, N);
rM = zeros(N,dim,rec_num); %record history positions
rM_cos = zeros(1,N,rec_num);
%% Main Loop
while FEs<MaxFEs
    calPositions = X;
    div_num = randperm(dim);
        %Divide the population and select the regions with more individuals based on the best
        for j=1:max(divide_num,1)
            th = best_pos(div_num(j));
            index = calPositions(:,div_num(j)) > th;
            if sum(index) < size(calPositions, 1)/2 %choose the side of the majority
                index = ~index;
            end
            calPositions = calPositions(index,:);
        end
    D = best_pos - calPositions; %Compute the distance from individuals to the best
    D_wind = sum(D, 1)/size(calPositions, 1); %Calculate the mean of all distances
 
    beta = size(calPositions, 1) / N;
    gama = 1/sqrt(1-power(beta,2));
    step = w * (rand(size(D_wind))-0.5) * (1-FEs/MaxFEs);
    step2 = 0.1*w*(rand(size(D_wind))-0.5)* (1-FEs/MaxFEs)*(1+1/2*(1+tanh(beta/gama))*(1-FEs/MaxFEs));
    step3 = 0.1*(rand()-0.5) * (1-FEs/MaxFEs);
    act =actCal(1 ./ 1 + (0.5 - 10*(rand(size(D_wind)))));
    
    if rec == 1 %record the first generation of positions
        rM(:,:,rec) = X;
        rM_cos(1,:,rec) = AllFitness;
        rec = rec + 1;
    end
  
    for i=1:N    
        newX(i,:) = X(i,:);
        %Spore dispersal search
        %Update M using Eq.(6)
         if rand()>d1
            newX(i,:) = newX(i,:) + step .* D_wind; %%1
         else
            newX(i,:) = newX(i,:) + step2 .* D_wind; %%2
        end

        if rand() < 0.8
            % Dual propagation search
            %Update M using Eq.(11)
            if rand() > 0.5
                newX(i,div_num(1)) = best_pos(div_num(1)) + step3 * D_wind(div_num(1)); %%3
            else
                newX(i,:) = (1-act) .* newX(i,:)+act .* best_pos; %%4
            end
        end
              
        %Boundary absorption
        Flag4ub=newX(i,:)>ub;
        Flag4lb=newX(i,:)<lb;
        newX(i,:)=(newX(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
 
        newAllFitness(i)=fobj(newX(i,:));
        FEs=FEs+1;
        
        % Cryptobiosis mechanism
        rM(i,:,rec) = newX(i,:);
        rM_cos(1,i,rec) = newAllFitness(i);
        
        
        if newAllFitness(i)<bestFitness 
           best_pos=newX(i,:);
           bestFitness=newAllFitness(i);
        end    
    end %End for
    
    rec = rec + 1;
    % Cryptobiosis mechanism
    if rec > rec_num || FEs>=MaxFEs
        [lcost,Iindex] = min(rM_cos, [] ,3);
        for i=1:N
            X(i,:) = rM(i,:,Iindex(i));
        end
        AllFitness = lcost;
        rec = 1;
    end
    
    Convergence_curve(it)=bestFitness;
    it=it+1;
end
end

function [act] = actCal(X)
    act = X;
    act(act>=0.5) = 1;
    act(act<0.5) = 0;
end

