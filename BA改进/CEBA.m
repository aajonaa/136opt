% Main programs starts here
function [best,Convergence_curve]=CEBA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
n=SearchAgents_no;  
A=0.9*ones(1,n);      % Loudness  (constant or decreasing)
r=0.5*ones(1,n);      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
% You should change these values if necessary
Qmin=0;         % Frequency minimum
Qmax=2;         % Frequency maximum
FEs=0;
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
  %Fitness(i)=fobj(Sol(i,:));
end
%%
%%混沌初始
x = rand();
while(x==0.25 || x==0.5 || x==0.75)
    x = rand();
end
ch(1) = x;
for ant=1:(n-1)
    ch(ant+1)=4*ch(ant)*(1-ch(ant));
    PCh = ch(ant)*Sol; 
    PHe = [Sol;PCh];
    count=size(PHe,1);
    FitnessHe1=[];
    for i=1:count
        Flag4ub=PHe(i,:)>ub;%限制边界
        Flag4lb=PHe(i,:)<lb;
        PHe(i,:)=(PHe(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;    
        PHeLin=fobj(PHe(i,:));
        FEs = FEs +1;
        FitnessHe1 = [FitnessHe1 PHeLin];
    end
%     FitnessHe1=calculateFitness(ObjValHe1);
    [FitnessHe2,index] = sort(FitnessHe1);
    Fitness = FitnessHe2(1:n);
    Sol = PHe(index,:);
    Sol = Sol(1:n,:);
end
%%
% Find the initial best solution
[fmin,I]=min(Fitness);
best=Sol(I,:);

Convergence_curve=[];
t=1;
GEN = 1;
% Main loop
while  FEs < MaxFEs
% for t=1:N_gen, 
% Loop over all bats/solutions\

        for i=1:n,
          Q=Qmin*ones(1,d)+(Qmin-Qmax)*rand(1,d);
          v(i,:)=0.2*v(i,:)+(Sol(i,:)-best).*Q;
          S(i,:)=Sol(i,:)+v(i,:);
          % Pulse rate
          if rand>r(i)
          % The factor 0.001 limits the step sizes of random walks 
              S(i,:)=best+randn(1,d);
          end

     % Evaluate new solutions
           Flag4ub=S(i,:)>ub;%限制边界
           Flag4lb=S(i,:)<lb;
           S(i,:)=(S(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            if FEs<MaxFEs
                FEs=FEs+1;
                Fnew=fobj(S(i,:));
               % Update if the solution improves, or not too loud
               if (Fnew<=Fitness(i)) && (rand<A(i))
                    Sol(i,:)=S(i,:);
                    Fitness(i)=Fnew;
                    A(i) = 0.9*A(i);
                    r(i) = r(i)*(1-exp((-0.9)*t));
                    if A(i) < 0.5
                        A(i) = 0.5;
                    end
                    if r(i) < 0.3
                       r(i)= 0.3;
                    end
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
        %%
        if GEN ==100
            %混沌初始
           x = rand();
           while(x==0.25 || x==0.5 || x==0.75)
               x = rand();
           end
           ch(1) = x;
            for ant=1:(n-1)
                ch(ant+1)=4*ch(ant)*(1-ch(ant));
                PCh = ch(ant)*Sol;
                PHe = PCh;
                count=size(PHe,1);
                FitnessHe1=[];
                for i=1:count
                    Flag4ub=PHe(i,:)>ub;%限制边界
                    Flag4lb=PHe(i,:)<lb;
                    PHe(i,:)=(PHe(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                    PHeLin=fobj(PHe(i,:));
                    FEs = FEs +1;
                    FitnessHe1 = [FitnessHe1 PHeLin];
                end
                PHe = [Sol;PCh];
                FitnessHe1 = [Fitness FitnessHe1];
                %     FitnessHe1=calculateFitness(ObjValHe1);
                [FitnessHe2,index] = sort(FitnessHe1);
                Fitness = FitnessHe2(1:n);
                Sol = PHe(index,:);
                Sol = Sol(1:n,:);
                if min(Fitness)<fmin
                    [fmin,I]=min(Fitness);
                    best=Sol(I,:);
                end
%                  %%
%                  best1 = Lb + Ub - best;
%                    FEs = FEs +1;
%                  if fobj(best1) < fmin
%                      fmin = fobj(best1);
%                      best = best1;
%                  end
            end
            
            GEN = 1;
        end
       
%%
        Convergence_curve(t)=fmin;
        t=t+1;
        GEN = GEN+1;
end

end

