function [Positionbest,Convergence_curve] = FA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
%  parameters [n N_iteration alpha betamin gamma]
% Simple bounds/limits for d-dimensional problems
Lb=lb.*ones(1,dim);
Ub=ub.*ones(1,dim);

u0=Lb+(Ub-Lb).*rand(1,dim);
% MaxGeneration=number of pseudo time steps
% ------------------------------------------------
% alpha=0.25;      % Randomness 0--1 (highly random)
% betamn=0.20;     % minimum value of beta
% gamma=1;         % Absorption coefficient
% ------------------------------------------------
n=SearchAgents_no;
alpha=0.5;
betamin=0.2;
gamma=1;

Lightbest=inf;
Convergence_curve=[];
% Calcualte dimension
dim=length(u0);
% Initial values of an array
zn=zeros(n,1);
[ns,Lightn]=init_ffa(n,dim,Lb,Ub,u0);

FEs=0;
t=1;
% Main loop
while  FEs < MaxFEs

% This line of reducing alpha is optional
 alpha=alpha_new(alpha,MaxFEs);

% Evaluate new solutions (for all n fireflies)
for i=1:n,
   if FEs<MaxFEs
       FEs=FEs+1;
       zn(i)=fobj(ns(i,:));
       Lightn(i)=zn(i);
   else
       break;
   end
end

% Ranking fireflies by their light intensity/objectives
% [Lightn,Index]=sort(zn,'descend');%------for max-------
[Lightn,Index]=sort(zn);
 ns=ns(Index,:);

%% Find the current best
nso=ns; Lighto=Lightn;
if Lightn(1)<Lightbest
    Lightbest=Lightn(1);%for min
    Positionbest=ns(1);
end
    Convergence_curve(t)=Lightbest;
    t=t+1;
    Positionbest=ns(1);

% Move all fireflies to the better locations
[ns]=ffa_move(n,dim,ns,Lightn,nso,Lighto,alpha,betamin,gamma,Lb,Ub);

end   %%%%% end of iterations

end

% -------------------------------------------------------
% ----- All the subfunctions are listed here ------------
% The initial locations of n fireflies
function [ns,Lightn]=init_ffa(n,d,Lb,Ub,u0)
  % if there are bounds/limits,
if ~isempty(Lb),
   for i=1:n,
   ns(i,:)=Lb+(Ub-Lb).*rand(1,d);
   end
else
   % generate solutions around the random guess
   for i=1:n,
   ns(i,:)=u0+randn(1,d);
   end
end

% initial value before function evaluations
Lightn=zeros(n,1);
end
% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,alpha,betamin,gamma,Lb,Ub)
% Scaling of the system
scale=abs(Ub-Lb);

% Updating fireflies
for i=1:n,
% The attractiveness parameter beta=exp(-gamma*r)
   for j=1:n,
      r=sqrt(sum((ns(i,:)-ns(j,:)).^2));
      % Update moves
        if Lightn(i)>Lighto(j), % Brighter and more attractive
            for dim=1:d
               beta0=1; beta=(beta0-betamin)*exp(-gamma*r.^2)+betamin;
               tmpf=alpha*(rand-0.5)*scale(1);
               ns(i,dim)=ns(i,dim)*(1-beta)+nso(j,dim)*beta+tmpf;
            end
        end
   end % end for j

end % end for i

% Check if the updated solutions/locations are within limits
[ns]=findlimits(n,ns,Lb,Ub);
end

% This function is optional, as it is not in the original FA
% The idea to reduce randomness is to increase the convergence,
% however, if you reduce randomness too quickly, then premature
% convergence can occur. So use with care.
function alpha=alpha_new(alpha,NGen)
% alpha_n=alpha_0(1-delta)^NGen=10^(-4);
% alpha_0=0.9
delta=1-(10^(-4)/0.9)^(1/NGen);
alpha=(1-delta)*alpha;
end
% Make sure the fireflies are within the bounds/limits
function [ns]=findlimits(n,ns,Lb,Ub)
for i=1:n,
     % Apply the lower bound
  ns_tmp=ns(i,:);
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);

  % Apply the upper bounds
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move
  ns(i,:)=ns_tmp;
end
end
%% ==== End of Firefly Algorithm implementation ======

