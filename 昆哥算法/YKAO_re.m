function [Leader_pos,Convergence_curve]=AO(N,Max_FEs,lb,ub,dim,fobj)

    FEs=0;
    it=1;
    Fitnorm=zeros(1,N);
    Convergence_curve=[];
    GV = zeros(N, dim);        

    X=initialization(N,dim,ub,lb);

    for i=1:N
        GV(i, :) = X(i, :) ./ (ub - lb);
        AllFitness(i)=fobj(X(i,:));
        FEs=FEs+1;
    end

    [fmin,x]=min(AllFitness);

    newX=zeros(N,dim);
    best=X(x,:);
    bestFitness=fmin;

    while FEs<=Max_FEs
        
        K= 1-((FEs)^(1/6)/(Max_FEs)^(1/6));

        E =1*exp(-4*(FEs/Max_FEs));
        lamda_t = 0.1 + (0.518 * ((1-(FEs/Max_FEs)^0.5))); 
    
        for i=1: N
            Fitnorm(i)= (AllFitness(i)-min(AllFitness))/(max(AllFitness)-min(AllFitness));
            for j=1:dim

                if rand<K 
                    if rand<0.5
                        newX(i,j) = X(i,j)+E.*X(i,j)*(-1)^FEs; 
                    else

                        newX(i,j) = X(i,j)+E.*best(j)*(-1)^FEs; %% 2

                    end
                else

                    newX(i,j)=X(i,j);

                end

                omega_it=(rand/2)+0.1;
                while true, r1 = round(N * rand + 0.5); if r1 ~= i && r1 ~= x, break, end, end
                while true, r2 = round(N * rand + 0.5); if r2 ~= i && r2 ~= x && r2 ~= r1, break, end, end
                newX(i,j) = best(j) + ((X(r2,j) - X(i,j)) * lamda_t) + ((X(r1,j) - X(i,j)) * omega_it);
                
            end


     
            ii = i + 1;
            if i == N
                ii = 1;
            end
            beta_1 = 1 + (rand / 2); 
            if  AllFitness(i) < beta_1 * bestFitness
                newX(i, :) = X(i, :) + abs(randn(1, dim)) .* (X(ii, :) - X(i, :)) + randn(1, dim) .* GV(i, :);
            end
            GV(i, :) = GV(i, :) .* ((rand ^ 2) * randn(1, dim));

            
   
            newX(i,:)=Transborder_reset(newX(i,:),ub,lb,dim,best); 
           
            tFitness=fobj(newX(i,:));
            FEs=FEs+1;
          
            if tFitness<AllFitness(i)
                X(i,:)= newX(i,:);
                AllFitness(i)=tFitness;
            end
        end
        [fmin,x]=min(AllFitness);
        if fmin<bestFitness
            best=X(x,:);
            bestFitness=fmin;
        end


        [~, SortOrder] = sort(AllFitness);
        
    
        SortOrder = SortOrder(1:N); 

        X = X(SortOrder, :);
        AllFitness = AllFitness(SortOrder);
        GV = GV(SortOrder, :);

       
        Convergence_curve(it)=bestFitness;
        Leader_pos=best;
        bestFitness = min(AllFitness);
        it=it+1;
    end
end

function z=Transborder_reset(z,ub,lb,dim,best)
    for j=1:dim
        if z(j)>ub || z(j)<lb
            
            z(j)=best(j);
            
        end
    end
end