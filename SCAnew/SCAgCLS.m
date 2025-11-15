function [Destination_position,Convergence_curve]=SCAgCLS(N,Max_iteration,lb,ub,dim,fobj)
tic
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=[];
Objective_values = zeros(1,size(X,1));
t=1;
 FEs=0;
% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
     FEs=FEs+1;
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end
Convergence_curve(t)=Destination_fitness;
%Main loop
t=t+1; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while  FEs<=Max_iteration
    
    % Eq. (3.4)
    a = 2;
    r1=a-t*((a)/Max_iteration);% r1 decreases linearly from a to 0
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
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
            
        end
    end
    
    for i=1:size(X,1)
         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate the objective values
        Objective_values(1,i)=fobj(X(i,:));
         f_X(i)=Objective_values(1,i);
         FEs=FEs+1;
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    
    %--------------------------------------------------%
       for i=1:N
        r=normrnd(0, 1, N,dim);
        t1=randperm(dim);
        cauchy=r(t1(1),:);
        X_meaution(i,:)=X(i,:)+rand(1,dim).*cauchy;
        Flag4ub=X_meaution(i,:)>ub;
        Flag4lb=X_meaution(i,:)<lb;
        X_meaution(i,:)=(X_meaution(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
         f_X_meaution(i)=fobj( X_meaution(i,:));
          FEs=FEs+1;
          if f_X_meaution(i) <  Destination_fitness
          Destination_position = X_meaution(i,:);
          Destination_fitness=f_X_meaution(i);
          end
      end 
         pos_combine = [X;X_meaution];
         total_meaution=[f_X,f_X_meaution];
           [~,q] = sort(total_meaution);
            for i = 1:N
           X(i,:) = pos_combine(q(i),:);
            end
    
  
   %--------------------------------------------------%
   [ Destination_position,Destination_fitness,FEs ] = CLS(N,Max_iteration,FEs,lb,ub,Destination_position,Destination_fitness,fobj );


    Convergence_curve(t)=Destination_fitness;
    
    % Increase the iteration counter
    t=t+1;
end
 toc
end
function r=C(x,y,N,dim)
    
Z=rand(N,dim);
r=x+y.*tan(pi.*(Z-0.5));

end
