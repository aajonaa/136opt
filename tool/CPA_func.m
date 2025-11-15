function [X,fes,count,conv] = CPA_func(N,Max_iteration,X,Target,Target_fitness,lb,ub,dim,fes,fobj)
count=0;
conv=[];

     a = exp(9-18*fes/Max_iteration);%公式5更新a   
    S0=a*(1-fes/Max_iteration); % r1 decreases linearly from a to 0

    
    %Communication and Collaboration
    for j=1 : dim
        EG=X(:,j);     %Selecting a dimensional population
        EG=sort(EG);
        i=round(rand()*(size(X,1))+0.5);  %Random selection of an independent individual
        if rand()< fes/Max_iteration
             X(i,:)=Target;
        end
            r=rand();
            X(i,j)=r*X(i,j)+(1-r)*(EG(1)+EG(2))/2; %Improve its position with two best individuals
    end
    

    for i=1:N % in i-th solution
        
        S=2*S0*rand-S0; %A represents the strength of the prey, which decreases with the number of iterations
        
        l=rand();
        
        if  abs(S) < 2*a/3 
            if  rand() > 0.5
                %Disperse food
                X(i,:)=Target-S*(rand(1,dim)*(ub(1)-lb(1))+lb(1));
            else
                %Encircle food
                Dp=abs(Target-X(i,:));
                X(i,:)=Target-2*S*Dp*exp(l)*tan(l*pi/4);
            end
        else
            
            distance=4*rand()-2;   %Represents distance from prey
            
            if abs(distance) < 1  %% [-2) -1...1 (2]  公式14
                %Supporting closest individual
                for j= 1:dim
                    P(j,:)=rand(1,dim)*X(i,j);
                    p_value(j,:)=fobj(P(j,:));
                    fes=fes+1;
                    count = count + 1;
                    conv(1,count) = Target_fitness;
                end
                [p_value,index]=sort(p_value);
                X(i,:)=P(index(1),:);
            else                                       %公式14
                %Research food
                X_rand=rand(1,dim).*(ub(1)-lb(1))+lb(1);
                D=abs(2*rand()*X_rand-X(i,:));
                %%
                X(i,:)=X_rand-S*D;
            end
        end
    end

end