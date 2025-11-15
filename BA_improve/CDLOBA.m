% Main programs starts here
function [best,Convergence_curve]=CDLOBA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
n=SearchAgents_no;  
A=rand(1,n)+ones(1,n);      % Loudness  (constant or decreasing)
r=rand(1,n);      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
% You should change these values if necessary
Qmin=0;         % Frequency minimum
Qmax=2;         % Frequency maximum

% Dimension of the search variables
d=dim;           % Number of dimensions 
% Lower limit/bounds/ a vector
Lb=lb.*ones(1,d);
% Upper limit/bounds/ a vector
Ub=ub.*ones(1,d);
% Initializing arrays
Q=zeros(n,1);   % Frequency
v=zeros(n,d);   % Velocities
% Initialize the population/solutions
for i=1:n,
  Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
  Fitness(i)=fobj(Sol(i,:));
end
S=Sol;
% Find the initial best solution
[fmin,I]=min(Fitness);
best=Sol(I,:);

Convergence_curve=[];
FEs=0;
t=1;
% Main loop
while  FEs < MaxFEs
    % for t=1:N_gen,
    % Loop over all bats/solutions
    FitnessHe1=[];
    for i = 1:n
        for j = 1:d
            S_opp(i,j)=(1/2)*(Lb(j)+Ub(j))-S(i,j);
        end
        finess_S_opp(i)=fobj(S_opp(i,:));
        finess_S(i)=fobj(S(i,:));
        FEs=FEs+1;
        FitnessHe1 = [FitnessHe1 finess_S_opp(i)];
    end
    FitnessHe1=[FitnessHe1 finess_S];
    PHe=[S_opp;S];
    [FitnessHe2,index] = sort(FitnessHe1);
    fitness = FitnessHe2(1:n);
    S1 = PHe(index,:);
    S = S1(1:n,:);
    if fitness(1)<fmin
        fmin=fitness(1);
        best=S(1,:);
    end
    for i=1:n,
        Q=Qmin+(Qmin-Qmax)*rand(1,d);
        v(i,:)= v(i,:)+(Sol(i,:)-best).*Q;
        S(i,:)=Sol(i,:)+v(i,:);
        % Apply simple bounds/limits
        Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
        % Pulse rate
        if rand>r(i)
            % The factor 0.001 limits the step sizes of random walks
            S(i,:)=best+A(i)*randn(1,d);
        end
        
        % Evaluate new solutions
        S(i,:)=simplebounds(S(i,:),Lb,Ub);
        if FEs<MaxFEs
            FEs=FEs+1;
            Fnew=fobj(S(i,:));
            % Update if the solution improves, or not too loud
            if (Fnew<=Fitness(i)) && (rand<A(i))
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
                A(i) = 0.9*A(i);
                r(i) = r(i)*(1-exp((-0.9)*t));
            end
            
        else
            break;
        end
        
        % Update the current best solution
        if Fnew<=fmin,
            best=S(i,:);
            fmin=Fnew;
        end
    end
    Convergence_curve(t)=fmin;
    t=t+1;
end

end
% Application of simple limits/bounds
function s=simplebounds(s,Lb,Ub)
% Apply the lower bound vector
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);

% Apply the upper bound vector
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
% Update this new move
s=ns_tmp;
end
