% Qu, C., et al. (2018). "A Modified Sine-Cosine Algorithm Based on Neighborhood Search and Greedy Levy Mutation." Computational Intelligence and Neuroscience 2018.

function [Destination_position, Convergence_curve]=MSCA2(N,MaxFEs,lb,ub,dim,fobj)
%A Modified Sine-Cosine Algorithm Based on Neighborhood Search and Greedy Levy Mutation
fes = 0;
%Set the initial parameters
wMax = 0.9;
wMin = 0.4;
epsilon = 30;
lambda = 0.01;
a = 2;

r = zeros(1,dim);
rMax = zeros(1,dim);
theta = zeros(1,dim);
beta = 1.5;
deltaU = ((gamma(1+beta)*sin(pi*beta/2)) / (gamma((1+beta)/2)*beta*(2^((beta-1)/2)))) ^ (1/beta);
deltaV = 1;
gs_best = zeros(1,dim);

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Objective_values = zeros(1,size(X,1));

% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
    fes = fes + 1;
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
end

%Main loop
t=1;
while fes<=MaxFEs
    
    w = wMax - (wMax - wMin)*fes/MaxFEs;
    
    r1 = a*exp(fes/MaxFEs);
    
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=(2 - fes/MaxFEs)*rand();
            r4=rand();
            
            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,j)= w*X(i,j)+r1*sin(r2)*abs(r3*Destination_position(j)*(1 + lambda*unifrnd(-1,1))-X(i,j));
            else
                % Eq. (3.2)
                X(i,j)= w*X(i,j)+r1*cos(r2)*abs(r3*Destination_position(j)*(1 + lambda*unifrnd(-1,1))-X(i,j));
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
        fes = fes + 1;
        
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    
    % optimal individual based on greedy levy variation
    for j = 1 : dim
        gs_best = Destination_position;
        r(j) = abs(Destination_position(j) - (1/N)*sum(X(:,j)));
        rMax(j) = max(X(:,j)) - min(X(:,j));
        theta(j) = exp((-epsilon*fes/MaxFEs)*(1 - r(j)/rMax(j)));
        levy = normrnd(0,deltaU^2)/(abs(normrnd(0,deltaV^2))^(1/beta));
        gs_best(j) = gs_best(j) + theta(j)*levy*Destination_position(j);
        gs_best_fit = fobj(gs_best);
        fes = fes + 1;
        if gs_best_fit < Destination_fitness
            Destination_position(j) = gs_best(j);
            Destination_fitness = gs_best_fit;
        end
    end
    
    Convergence_curve(t)=Destination_fitness;
    
    % Increase the iteration counter
    t=t+1;
end
end
