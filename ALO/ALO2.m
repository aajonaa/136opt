% ALO modified by Jona 2024-3-12.

function [Elite_antlion_position,Convergence_curve]=ALO(N,Max_FEs,lb,ub,dim,fobj)
tic
% disp('ALO');
% Initialize the positions of antlions and ants
antlion_position=initialization(N,dim,ub,lb);
ant_position=initialization(N,dim,ub,lb);

FEs = 0;
iter = 1;

% Initialize variables to save the position of elite, sorted antlions, 
% convergence curve, antlions fitness, and ants fitness
Sorted_antlions=zeros(N,dim);
Elite_antlion_position=zeros(1,dim);
Elite_antlion_fitness=inf;
Convergence_curve=[];
antlions_fitness=zeros(1,N);
ants_fitness=zeros(1,N);

% Calculate the fitness of initial antlions and sort them
for i=1:size(antlion_position,1)
    antlions_fitness(1,i)=fobj(antlion_position(i,:)); 
    FEs = FEs + 1;
end

[sorted_antlion_fitness,sorted_indexes]=sort(antlions_fitness);
    
for newindex=1:N
     Sorted_antlions(newindex,:)=antlion_position(sorted_indexes(newindex),:);
end
    
Elite_antlion_position=Sorted_antlions(1,:);
Elite_antlion_fitness=sorted_antlion_fitness(1);

% Main loop start from the second iteration since the first iteration 
% was dedicated to calculating the fitness of antlions

while FEs < Max_FEs
    disp(FEs);
%     This for loop simulate random walks
    for i=1:size(ant_position,1)
        % Select an antlion based on their fitness (the better anlion the higher chance of catching ant)
        Rolette_index=RouletteWheelSelection(1./sorted_antlion_fitness)
        if Rolette_index==-1  
            Rolette_index=1;
        end
      
        % RA is the random walk around the selected antlion by rolette wheel
        antlion = Sorted_antlions(Rolette_index,:)
        if size(lb,1) ==1 && size(lb,2)==1 %Check if the bounds are scalar
            lb=ones(1,dim)*lb;
            ub=ones(1,dim)*ub;
        end

        if size(lb,1) > size(lb,2) %Check if boundary vectors are horizontal or vertical
            lb=lb';
            ub=ub';
        end

        I=1; % I is the ratio in Equations (2.10) and (2.11)

        if iter>Max_FEs/dim/10
            I=1+100*(iter/Max_FEs/dim)
        end

        if iter>Max_FEs/dim/2
            I=1+1000*(iter/Max_FEs/dim)
        end

        if iter>Max_FEs/dim*(3/4)
            I=1+10000*(iter/Max_FEs/dim)
        end

        if iter>Max_FEs/dim*(0.9)
            I=1+100000*(iter/Max_FEs/dim)
        end

        if iter>Max_FEs/dim*(0.95)
            I=1+1000000*(iter/Max_FEs/dim)
        end


        % Dicrease boundaries to converge towards antlion
        lb=lb/(I) % Equation (2.10) in the paper
        ub=ub/(I) % Equation (2.11) in the paper

        % Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
        if rand<0.5
            lb=lb+antlion % Equation (2.8) in the paper
        else
            lb=-lb+antlion
        end

        if rand>=0.5
            ub=ub+antlion % Equation (2.9) in the paper
        else
            ub=-ub+antlion
        end

        % This function creates n random walks and normalize accroding to lb and ub
        % vectors 
        for i=1:dim
            X = [0 cumsum(2* (rand(Max_FEs/dim,1)>0.5) -1)']; % Equation (2.1) in the paper
            %[a b]--->[c d]
            a=min(X)
            b=max(X)
            c=lb(i)
            d=ub(i)     
            X_norm=((X-a).*(d-c))./(b-a)+c; % Equation (2.7) in the paper
            RWs(:,i)=X_norm;
            disp(size(RWs))
        end
%         RA=Random_walk_around_antlion(dim,Max_FEs/dim,lb,ub, Sorted_antlions(Rolette_index,:),iter);
        RA = RWs(1, :)
        
        % RE is the random walk around the elite (best antlion so far)
        antlion = Sorted_antlions(1,:)
        if size(lb,1) ==1 && size(lb,2)==1 %Check if the bounds are scalar
            lb=ones(1,dim)*lb;
            ub=ones(1,dim)*ub;
        end

        if size(lb,1) > size(lb,2) %Check if boundary vectors are horizontal or vertical
            lb=lb';
            ub=ub';
        end

        I=1; % I is the ratio in Equations (2.10) and (2.11)

        if iter>Max_FEs/dim/10
            I=1+100*(iter/Max_FEs/dim)
        end

        if iter>Max_FEs/dim/2
            I=1+1000*(iter/Max_FEs/dim)
        end

        if iter>Max_FEs/dim*(3/4)
            I=1+10000*(iter/Max_FEs/dim)
        end

        if iter>Max_FEs/dim*(0.9)
            I=1+100000*(iter/Max_FEs/dim)
        end

        if iter>Max_FEs/dim*(0.95)
            I=1+1000000*(iter/Max_FEs/dim)
        end


        % Dicrease boundaries to converge towards antlion
        lb=lb/(I) % Equation (2.10) in the paper
        ub=ub/(I) % Equation (2.11) in the paper

        % Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
        if rand<0.5
            lb=lb+antlion % Equation (2.8) in the paper
        else
            lb=-lb+antlion
        end

        if rand>=0.5
            ub=ub+antlion % Equation (2.9) in the paper
        else
            ub=-ub+antlion
        end

        % This function creates n random walks and normalize accroding to lb and ub
        % vectors 
        for i=1:dim
            X = [0 cumsum(2* (rand(Max_FEs/dim,1)>0.5) -1)']; % Equation (2.1) in the paper
            %[a b]--->[c d]
            a=min(X)
            b=max(X)
            c=lb(i)
            d=ub(i)  
            X_norm=((X-a).*(d-c))./(b-a)+c; % Equation (2.7) in the paper
            RWs(:,i)=X_norm;
            disp(size(RWs));
        end
%         RE=Random_walk_around_antlion(dim,Max_FEs/dim,lb,ub, Elite_antlion_position(1,:),iter);
        RE = RWs(1, :)

%         ant_position(i,:)= (RA(iter,:)+RE(iter,:))/2; % Equation (2.13) in the paper  
        ant_position(i, :) = (RA + RE) / 2
    end
    
    for i=1:size(ant_position,1)  
        
        % Boundar checking (bring back the antlions of ants inside search
        % space if they go beyoud the boundaries
        Flag4ub=ant_position(i,:)>ub;
        Flag4lb=ant_position(i,:)<lb;
        ant_position(i,:)=(ant_position(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        
        ants_fitness(1,i)=fobj(ant_position(i,:));  
        FEs = FEs + 1;
       
    end
    
    % Update antlion positions and fitnesses based of the ants (if an ant 
    % becomes fitter than an antlion we assume it was cought by the antlion  
    % and the antlion update goes to its position to build the trap)
    double_population=[Sorted_antlions;ant_position];
    double_fitness=[sorted_antlion_fitness ants_fitness];
        
    [double_fitness_sorted, I]=sort(double_fitness);
    double_sorted_population=double_population(I,:);
        
    antlions_fitness=double_fitness_sorted(1:N);
    Sorted_antlions=double_sorted_population(1:N,:);
        
    % Update the position of elite if any antlinons becomes fitter than it
    if antlions_fitness(1)<Elite_antlion_fitness 
        Elite_antlion_position=Sorted_antlions(1,:);
        Elite_antlion_fitness=antlions_fitness(1);
    end
      
    % Keep the elite in the population
    Sorted_antlions(1,:)=Elite_antlion_position;
    antlions_fitness(1)=Elite_antlion_fitness;
  
    % Update the convergence curve
    Convergence_curve(iter)=Elite_antlion_fitness;
    iter=iter+1; 
end
toc
end






