% 混沌局部搜索 CLS + SCA
function [Convergence_curve]=CLSSCA(N,Max_iteration,lb,ub,dim,fobj,funcNum)
fobj = str2func(['cec14_func_', num2str(funcNum)]);
tic
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
    
    All_objective_values(1,i)=Objective_values(1,i);
end

%Main loop
t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while t<=Max_iteration
    
    % Eq. (3.4)
    a = 2;
    %     Max_iteration = Max_iteration;
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
        %% 混沌局部搜索
        setCan = (Max_iteration-t+1)/Max_iteration;
        %         setCan = 1- power(abs((t-1)/t),1000);
        x = rand();
        while(~(x~=0.25 && x~=0.5 && x~=0.75 && x~=1))
            x=rand();
        end
        ch(1) = x;
        for ant=1:N
            ch(ant+1)=4*ch(ant)*(1-ch(ant));
            CH(ant) = lb+ch(ant)*(ub-lb);    %ub大
            V = (1-setCan)*Destination_position+setCan*CH(ant);
            % ObjValV=feval(objfun,V);         %计算函数值

            %% 边界控制
            Flag4ub=V>ub';
            Flag4lb=V<lb';
            V=(V.*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
            %% 边界控制结束

            FitnessV=fobj(V);%计算适应度值
            %    if (ObjValV<GlobalMin)
            if (FitnessV<Destination_fitness)
                Destination_fitness = FitnessV;
                Destination_position = V;
                %             break;
            end
        end
        %% 混沌最优值结束
        
        %% 高斯变异开始
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Moth_pos_m_gaus=X(i,:)*(1+randn(1));
        Moth_fitness_m_gaus=fobj(Moth_pos_m_gaus);
        Moth_fitness_s=fobj(X(i,:));

        Moth_fitness_comb=[Moth_fitness_m_gaus,Moth_fitness_s];
        [~,m]=min(Moth_fitness_comb);
        if m==1
            X(i,:)=Moth_pos_m_gaus;
        else
            X(i,:)=X(i,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 高斯变异结束
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
toc
Convergence_curve(1)=Convergence_curve(2);

end

function o=Levy(d)
beta=1.5;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
% Eq. (3.9)
o=step;
end