% Main programs starts here
function [best,Convergence_curve]=RCBA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
n=SearchAgents_no;
% This frequency range determines the scalings
% You should change these values if necessary
Qmin=0;         % Frequency minimum
Qmax=2;         % Frequency maximum
A=rand(1,n);      % Loudness  (constant or decreasing)
%r=rand(1,n);      % Pulse rate (constant or decreasing)
r=ones(1,n)*0.5;
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
    u = randn();
    p = rand();
    for i=1:n,
        Q=Qmin+(Qmin-Qmax)*rand;
        v(i,:)= v(i,:)+(Sol(i,:)-best)*Q;
        S(i,:)=Sol(i,:)+v(i,:);
        % Apply simple bounds/limits
        Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
        % Pulse rate
        if rand>r(i)
            % The factor 0.001 limits the step sizes of random walks
            for j = 1:d
                if rand<=p
                    S(i,j)=best(1,j)+u;
                else
                    S(i,j)=best(j)+A(i)*randn;
                end
            end
        end
        
        % Evaluate new solutions
        S(i,:)=simplebounds(S(i,:),Lb,Ub);
        if FEs<MaxFEs
            FEs=FEs+1;
            Fnew=fobj(S(i,:));
            % Update if the solution improves, or not too loud
            if (Fnew<=Fitness(i))
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
            end
        else
            break;
        end
        
        % Update the current best solution
        if Fnew<=fmin,
            best=S(i,:);
            fmin=Fnew;
        end
        if A(i)<0.7
            A(i) = A(i)/0.7;
        else
            A(i)=(10*(1-A(i)))/3;
        end
        %r(i) = r(i)+0.2-(0.5/(2*pi))*sin(2*pi*r(i));
         r(i) = r(i)*(1-exp((-0.9)*t));
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
