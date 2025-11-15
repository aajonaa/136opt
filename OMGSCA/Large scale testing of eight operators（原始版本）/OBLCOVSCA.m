function [Convergence_curve]=OBLCOVSCA(N,Max_iteration,lb,ub,dim,fobj,funcNum)
fobj = str2func(['cec14_func_', num2str(funcNum)]);
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));

VPosition=zeros(N,dim);
VFitness=zeros(1,N);
RouNum=ceil(N/4);
RouBestPosition=zeros(RouNum,dim);
% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end

%Main loop
t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while t<=Max_iteration
    
    % Eq. (3.4)
    a = 2;
%     Max_iteration = Max_iteration;
    r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
      [~,Sort_Index]=sort(Objective_values);
    for i=1:RouNum
        RouBestPosition(i,:)=X(Sort_Index(i),:);
    end
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            
            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                VPosition(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
               VPosition(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
          VFitness(1,i)=fobj(VPosition(i,:));
        if VFitness(1,i)<Objective_values(1,i)
            Objective_values(1,i)=VFitness(1,i);
            X(i,:)=VPosition(i,:);
        end
        end
       
    end
    
%      [P] = OBL( X,N,dim,ub,lb,fobj,0);
     [P] = OBL( X,N,dim,ub,lb,fobj,0);
     X = P;
        for i=1:size(X,1)
        [ ZPosition ] = Covariance( RouBestPosition );
        if fobj(ZPosition)<Objective_values(1,i)
            Objective_values(1,i)=fobj(ZPosition);
            X(i,:)=ZPosition;
        end
        end
    
    for i=1:size(X,1)
         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate the objective values
        Objective_values(1,i)=fobj(X(i,:));
        
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    %%
%     if mod(t,5)==0
%         [ BestFitness,BestPosiiton ] = MRM( Destination_position,X,N,dim,fobj );
%         if BestFitness<FoodFitness
%             Destination_position=BestPosiiton;
%             X=BestFitness;
%         end
%     end
% 	[ best_pos, fmin ] = simplex( Destination_position, fobj);
% 	if fmin<Destination_fitness
% 		Destination_position=best_pos;                      
% 		Destination_fitness=fmin;	
% 	end
%    [~,BestIndex]=min(Destination_fitness);
%     [ k1,k2 ] = GetRan2( BestIndex,N );
%     [ CanPosition,CanFitness ] = QuaIn( X(k1,:),Destination_position,X(k2,:),fobj,dim);
%     if CanFitness<Destination_fitness
%         Destination_position=CanPosition;
%         Destination_fitness=CanFitness;
%     end
    %%
%        for i=1:size(SalpPositions,1) 
% %         for j=1:dim
% %             T1(j)=rand()*(ub(j)+lb(j)) - SalpPositions(i,j);
% %         end
%     %%T=Xk+rand()*(Food-Xk)
% %         [ k ] = GetRan( i,N );
% %         for j=1:dim
% %             T2(j)=SalpPositions(k,j)+rand()*(FoodPosition(j)-SalpPositions(k,j));
% %         end
% %         [ BestPosition,BestFitness ] = Con_OEDCanSolution( SalpPositions(i,:),T2,dim,Q,F,fobj );
%     %%T=Xk1+rand()*(Xk2-Xk3)
%        [ k1,k2,k3 ] = GetRan3( i,N );
%         for j=1:dim
%             T2(j)=X(k1,j)+rand()*(X(k2,j)-X(k3,j));
%         end
%         [ BestPosition,BestFitness ] = Con_OEDCanSolution( X(i,:),T2,dim,Q,F,fobj );
% %        [ BestPosition,BestFitness ] = Con_OEDCanSolution( SalpPositions(i,:),FoodPosition,dim,Q,F,fobj );
%         if BestFitness<fobj(X(i,:))
%             X(i,:)=BestPosition;
%         end
% %         if fobj(T1)<fobj(SalpPositions(i,:))
% %             SalpPositions(i,:)=T1;
% %         end
%        end
%    [ BestPosition,BestFitness] = ADN(Cr,dim,fobj)

    Convergence_curve(t)=Destination_fitness;
    
    % Increase the iteration counter
    t=t+1;
end
 Convergence_curve(1)=Convergence_curve(2);
end
