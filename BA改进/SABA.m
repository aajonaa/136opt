% Main programs starts here
function [Convergence_curve]=SABA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
n=SearchAgents_no;  
% A=rand;      % Loudness  (constant or decreasing)
% r=rand+1;      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
% You should change these values if necessary
Cw = 3;
fmax = 2.5;
% Dimension of the search variables
d=dim;           % Number of dimensions 
% Lower limit/bounds/ a vector
Lb=lb.*ones(1,d);
% Upper limit/bounds/ a vector
Ub=ub.*ones(1,d);
for i = d
    s1 = Ub(i)-Lb(i);
end
s = s1/d;
% Initializing arrays
v=zeros(n,d);   % Velocities
FEs=0;
% Initialize the population/solutions
for i=1:n,
  S(i,:)=Lb+(Ub-Lb).*rand(1,d);
  Fitness(i)=fobj(S(i,:));
  FEs = FEs +1;
end
% Find the initial best solution
[fmin,I]=min(Fitness);
best=S(I,:);
Favg = Fitness;Favg_Sol = S; Fbest = min(Fitness);Fbest_Sol = best;%记录个体最优值和全局最优值
Convergence_curve=[];
t=1;
% Main loop
while  FEs < MaxFEs
% for t=1:N_gen, 
% Loop over all bats/solutions
 k = FEs/MaxFEs;
 r1 = rand+0.5;r2 = rand+0.5;
 w = (FEs/MaxFEs)+0.1;
 u = (FEs/MaxFEs)*0.4+0.3;
 if w > 0.9
     w = 0.9;
 end    
        for i=1:n,
     
          f1 = 1-exp(-abs(Favg(i)-Fbest))+1.5*(1-k)+0.5;
          f2 = Cw-f1;
          A(i) = f1/fmax;R(i)=f2/fmax;
          v(i,:)=w*v(i,:)+f1*r1*(Favg_Sol(i,:)-S(i,:))+f2*r2* (Fbest_Sol-S(i,:));
          S(i,:)=S(i,:)+u*v(i,:);
          % Apply simple bounds/limits
          S(i,:)=simplebounds(S(i,:),Lb,Ub);
          % Pulse rate
          if rand*0.7<R(i)
          % The factor 0.001 limits the step sizes of random walks 
          if (k>0) && (k<=0.1)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*2;
          elseif (k>0.1) && (k<=0.2)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*1.5;
          elseif (k>0.2) && (k<=0.3)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*1;
          elseif (k>0.3) && (k<=0.4)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*0.5;
          elseif (k>0.4) && (k<=0.6)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*0.1;
          elseif (k>0.6) && (k<=0.7)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*(0.1^3);
          elseif (k>0.7) && (k<=0.8)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*(0.1^5);
          elseif (k>0.8) && (k<=0.9)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*(0.1^7);
         elseif (k>0.9) && (k<=1)
              S(i,:)=Fbest_Sol+A(i)*(2*rand(1,d)-ones(1,d))*s*(0.1^9);
          end
          elseif (rand>0.5) && (rand<A(i))
              S(i,:)=Lb+(Ub-Lb).*rand(1,d);
          end

     % Evaluate new solutions
           S(i,:)=simplebounds(S(i,:),Lb,Ub);
            if FEs<MaxFEs
                FEs=FEs+1;
                Fnew=fobj(S(i,:));
                if Fnew < Favg(i)
                    Favg_Sol(i,:) = S(i,:);
                    Favg(i) = Fnew;
                end    
%                % Update if the solution improves, or not too loud
%                if (rand>0.5) && (rand<A(i)) ,
%                     S(i,:)=Lb+(Ub-Lb).*rand(1,d);
%                     Fitness(i)=fobj(S(i,:));
%                end
            else
                break;
            end
            
          % Update the current best solution
          if Fnew<=fmin,
                best=S(i,:);
                fmin=Fnew;
          end
          if Fbest > fmin
              Fbest_Sol = best;
              Fbest = fmin;
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
