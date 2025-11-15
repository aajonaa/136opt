%ZQoblº¯Êý + OL +  SCA 
function [Convergence_curve]=OLSCAV2(N,Max_iteration,lb,ub,dim,fobj, funcNum)
fobj = str2func(['cec14_func_', num2str(funcNum)]);
tic
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
LB=ones(1,dim).*lb;
UB=ones(1,dim).*ub;
Destination_position=zeros(1,dim);
Destination_fitness=inf;
   SearchAgents_no = N;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));
V=zeros(N,dim);
VFitness=zeros(1,N);
Trial=zeros(1,N);
Limits=20;
Q=5;
F=6;
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
    Max_iteration = Max_iteration;
    r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
    
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
                V(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                V(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
            VFitness(1,i)=fobj(V(i,:));
            if VFitness<Objective_values(1,i)
                X(i,:)=V(i,:);
                Objective_values(1,i)=VFitness(1,i);
            else
                Trial(1,i)=Trial(1,i)+1;
            end
        end
        
        
%                  [ k1,k2,k3 ] = GetRan3( i,N );
%                 for j=1:dim
%                     T2(j)=X(k1,j)+rand()*(X(k2,j)-X(k3,j));
%                 end
%                 [ BestPosition,BestFitness ] = Con_OEDCanSolution( X(i,:),T2,dim,Q,F,fobj );
%         %        [ BestPosition,BestFitness ] = Con_OEDCanSolution( SalpPositions(i,:),FoodPosition,dim,Q,F,fobj );
%                 if BestFitness<fobj(X(i,:))
%                     X(i,:)=BestPosition;
%                 end
    end
    for i=1:size(X,1)
        if Trial(1,i)>=Limits
            [ BestPosition,BestFitness ] = Con_OEDCanSolution( X(i,:),Destination_position,dim,Q,F,fobj );
            if BestFitness<Objective_values(1,i)
                X(i,:)=BestPosition;
                Objective_values(1,i)=BestFitness;
            else
                for j=1:size(X,2)
                    X(i,j)=LB(j)+rand()*(UB(j)-LB(j));
                end
            end
            Trial(1,i)=0;
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
   
    
    Convergence_curve(t)=Destination_fitness;
    
    % Increase the iteration counter
    t=t+1;
end
Convergence_curve(1)=Convergence_curve(2);
toc
end
