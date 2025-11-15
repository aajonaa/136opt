%%  GB¹Ç¼Ü + SCA
function [Convergence_curve]=GBSCA(N,Max_iteration,lb,ub,dim,fobj,funcNum)
fobj = str2func(['cec14_func_', num2str(funcNum)]);
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));
V=zeros(N,dim);
Vfitness=zeros(1,N);
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
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
            
        end
    end
    
    
    %%¹Ç¼Ü±äÒì
    for i=1:N
        Objective_values(1,i)=fobj(X(i,:));
        %           CR = 1-exp(-abs(Objective_values(1,i)-Destination_fitness));%%±äÒì¿ØÖÆ±äÁ¿
        CR = 0.3;
        [ k1,k2 ] = GetRan2( i,N );
        k=rand();
        for j=1:dim
            if rand()<CR
                mu=(Destination_position(j)+X(i,j))/2;
                sigma=abs(Destination_position(j)-X(i,j));
                V(i,j)=normrnd(mu,sigma);
            else
                V(i,j)=Destination_position(j)+k*(X(k1,j)-X(k2,j));%%DE/best/1
            end
        end
        Vfitness(i)=fobj(V(i,:));
        if Vfitness(i)<Objective_values(1,i)
            X(i,:)=V(i,:);
        end
    end
    %%¹Ç¼Ü±äÒì½áÊø
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
end
