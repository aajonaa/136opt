
% Particle Swarm Optimization
function [bestPos,cg_curve]=PSO(N,MaxFEs,lb,ub,dim,fobj)

%PSO Infotmation

Vmax=6;
noP=N;
wMax=0.9;
wMin=0.2;
c1=2;
c2=2;

% Initializations
vel=zeros(noP,dim);
pBestScore=zeros(noP);
pBest=zeros(noP,dim);
gBest=zeros(1,dim);
cg_curve=[];

% Random initialization for agents.
pos=initialization(noP,dim,ub,lb); 

for i=1:noP
    pBestScore(i)=inf;
end

% Initialize gBestScore for a minimization problem
 gBestScore=inf;
     
it=1;    
FEs=1;
% Main loop
while  FEs < MaxFEs
    
    % Return back the particles that go beyond the boundaries of the search
    % space
     Flag4ub=pos(i,:)>ub;
     Flag4lb=pos(i,:)<lb;
     pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    
    for i=1:size(pos,1)     
        %Calculate objective function for each particle
        if FEs<MaxFEs
            FEs=FEs+1;    
            fitness=fobj(pos(i,:));

            if(pBestScore(i)>fitness)
                pBestScore(i)=fitness;
                pBest(i,:)=pos(i,:);
            end
            if(gBestScore>fitness)
                gBestScore=fitness;
                gBest=pos(i,:);
            end
        else
            break;
        end
    end

    %Update the W of PSO
%     w=wMax-l*((wMax-wMin)/iter);
    w=1;
    %Update the Velocity and Position of particles
    for i=1:size(pos,1)
        for j=1:size(pos,2)       
            vel(i,j)=w*vel(i,j)+c1*rand()*(pBest(i,j)-pos(i,j))+c2*rand()*(gBest(j)-pos(i,j));
            
            if(vel(i,j)>Vmax)
                vel(i,j)=Vmax;
            end
            if(vel(i,j)<-Vmax)
                vel(i,j)=-Vmax;
            end            
            pos(i,j)=pos(i,j)+vel(i,j);
        end
    end
    cg_curve(it)=gBestScore;
    it=it+1;
    bestPos=gBest;
end

end