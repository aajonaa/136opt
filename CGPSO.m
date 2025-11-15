
% Particle Swarm Optimization
function [bestPos,cg_curve]=CGPSO(N,MaxFEs,lb,ub,dim,fobj)
tic
%PSO Infotmation

Vmax=6;
noP=N;
wMax=0.9;
wMin=0.2;
c1=2;
c2=2;

% Initializations
iter=MaxFEs;
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
     
FEs=0;
it=1;
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
     %%%=================CLS=================================
    if FEs/MaxFEs<0.8
%         setCan = (MaxFEs-FEs+1)/MaxFEs;
        setCan = 1-power(abs((FEs-1)/FEs),1500);

        x = rand();
        if(x~=0.25&&x~=0.5&&x~=0.75)
            ch(1) = x;
        end
        for ant=1:(noP)
            ch(ant+1)=4*ch(ant)*(1-ch(ant));
            CH(ant,:) = lb+ch(ant)*(ub-lb);    %ub大
            V = (1-setCan)*gBest+setCan*CH(ant);
            if FEs<MaxFEs
                FEs=FEs+1;
                FitnessV=fobj(V);%计算适应度值
                if (FitnessV<gBestScore)
                    gBestScore = FitnessV;
                    gBest = V;
                    break;
                end
            else
                break;
            end
                
        end
    else 
        setCan = (MaxFEs-FEs+1)/MaxFEs;
%         
        ch=randn(1,noP);

        for ant=1:(noP)
            CH(ant,:) = lb+ch(ant)*(ub-lb);    %ub大
            V = (1-setCan)*gBest+setCan*CH(ant);
            if FEs<MaxFEs
                % FEs=FEs+1;
                FitnessV=fobj(V);%计算适应度值
                if (FitnessV<gBestScore)
                    gBestScore = FitnessV;
                    gBest = V;
                    break;
                end
            else
                break;
            end
        end
    end
 %%%=======================================================   
    cg_curve(it)=gBestScore;
    it=it+1;
    bestPos=gBest;
end
toc
end