% function [FoodFitness,FoodPosition,Convergence_curve]=SSA(N,Max_iter,lb,ub,dim,fobj)
function [FoodPosition,Convergence_curve]=SSA(N,MaxFEs,lb1,ub1,dim,fobj)

ub=ones(dim,1).*ub1';
lb=ones(dim,1).*lb1';


Convergence_curve = [];
FEs=0;
%Initialize the positions of salps
SalpPositions=initialization(N,dim,ub,lb);


FoodPosition=zeros(1,dim);
FoodFitness=inf;


%calculate the fitness of initial salps

for i=1:size(SalpPositions,1)
    SalpFitness(1,i)=fobj(SalpPositions(i,:));
    FEs=FEs+1;
end

[sorted_salps_fitness,sorted_indexes]=sort(SalpFitness);

for newindex=1:N
    Sorted_salps(newindex,:)=SalpPositions(sorted_indexes(newindex),:);
end

FoodPosition=Sorted_salps(1,:);
FoodFitness=sorted_salps_fitness(1);

%Main loop
l=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness of salps
while FEs<MaxFEs+2
    
    c1 = 2*exp(-(4*FEs/MaxFEs)^2); % Eq. (3.2) in the paper
    
    for i=1:size(SalpPositions,1)
        
        SalpPositions= SalpPositions';
        
        if i<=N/2
            for j=1:1:dim
                c2=rand();
                c3=rand();
                %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                if c3<0.5 
                    SalpPositions(j,i)=FoodPosition(j)+c1*((ub(j)-lb(j))*c2+lb(j));
                else
                    SalpPositions(j,i)=FoodPosition(j)-c1*((ub(j)-lb(j))*c2+lb(j));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        elseif i>N/2 && i<N+1
            point1=SalpPositions(:,i-1);
            point2=SalpPositions(:,i);
            
            SalpPositions(:,i)=(point2+point1)/2; % % Eq. (3.4) in the paper
        end
        
        SalpPositions= SalpPositions';
    end
        lb1=lb1';ub1=ub1';
    for i=1:size(SalpPositions,1)
        
        Tp=SalpPositions(i,:)>ub1';Tm=SalpPositions(i,:)<lb1';SalpPositions(i,:)=(SalpPositions(i,:).*(~(Tp+Tm)))+ub1'.*Tp+lb1'.*Tm;
        
        SalpFitness(1,i)=fobj(SalpPositions(i,:));
            FEs=FEs+1;
        if SalpFitness(1,i)<FoodFitness
            FoodPosition=SalpPositions(i,:);
            FoodFitness=SalpFitness(1,i);
            
        end
    end
    lb1=lb1';ub1=ub1';

    Convergence_curve(l-1)=FoodFitness;
    l = l + 1;
end
end

function o=Levy1(d)
beta=3/2;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);

% Eq. (3.9)
o=step;

end
%lamda 0.5/0.75/1.5
function o=Levy2(d)
beta=0.75;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);

% Eq. (3.9)
o=step;

end
function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,1); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end






