function [Destination_position,Convergence_curve]=m_SCA(N,Max_iteration,lb,ub,dim,fobj)
Max_iteration = 1000;
fes = 0;
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));

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
end
%Initialize the algorithm parameters ¨C N (size of population of solutions), T (maximum number of iterations as termination criteria), A , perturbation or jumping rateJ R and self-adaptation rate S R
JR = 0.1;
%Main loop
t=1; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while t <= Max_iteration
    
    display('m_SCA');
    display(t);
    
    % Eq. (3.4)
    a = 2;
    r1=a-a*(t/Max_iteration); % r1 decreases linearly from a to 0
    
    if (rand < JR)
        [P,fes] = OBL_For_MaxFEs( X,N,dim,ub,lb,fobj,0,fes );
        X = P;
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
    else
        % Update the position of solutions with respect to destination
        for i=1:size(X,1) % in i-th solution
            for j=1:size(X,2) % in j-th dimension
                
                % Update r2, r3, and r4 for Eq. (3.3)
                r2=(2*pi)*rand();
                r3=2*rand;
                r4=rand();
                SR=rand();

                % Eq. (3.3)
                if r4<0.5
                    % Eq. (3.1)
                    X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j))) + SR*(Destination_position(j) - X(i,j));
                else
                    % Eq. (3.2)
                    X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j))) + SR*(Destination_position(j) - X(i,j));
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

            % Update the destination if there is a better solution
            if Objective_values(1,i)<Destination_fitness
                Destination_position=X(i,:);
                Destination_fitness=Objective_values(1,i);
            end
        end
    end
    Convergence_curve(t)=Destination_fitness;
    
    % Increase the iteration counter
    t=t+1;
end
end
