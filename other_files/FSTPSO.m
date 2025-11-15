% Particle Swarm Optimization
function [bestPos,cg_curve]= FSTPSO(N,maxFES,lb,ub,dim,fobj)
% 基于模糊策略自调整的粒子群算法
% 来源：Fuzzy Self-Tuning PSO_ A settings-free algorithm for global optimization

%PSO Infotmation

Vmax=6;
Vmin=-6; 
% noP = floor(10 + 2 * sqrt(dim));
noP = N;

wMax=0.9;
wMin=0.2;
c1=2;
c2=2;
w = 0.9;

% Initializations
pBestScore=zeros(noP, 1);
pScoreL = zeros(noP, 1);
pBest=zeros(noP,dim);
pPosL = zeros(noP, dim);
gBest=zeros(1,dim);
fitness = zeros(noP, 1);
dmax = sqrt( dim * (ub - lb)^2 );
convergence = [];
% Random initialization for agents.
pos=initialization(noP,dim,ub,lb); 
vel=rand(noP,dim) * (Vmax - Vmin) + Vmin;
for i=1:noP
    pBestScore(i)=inf;
end
gBestScore=inf;
FES = 0;     

l=1;
while FES <= maxFES
    for i=1:size(pos,1)     
        %Calculate objective function for each particle
        pScoreL(i) = fitness(i);
        fitness(i)=fobj(pos(i,:));
        FES = FES+1;
        if(pBestScore(i)>fitness(i))
            pBestScore(i)=fitness(i);
            pBest(i,:)=pos(i,:);
        end
        if(gBestScore>fitness(i))
            gBestScore=fitness(i);
            gBest=pos(i,:);
        end
    end
    fdelta = max(fitness);
    
    %Update the Velocity and Position of particles
    for i=1:size(pos,1)
        if ( l >= 2 )
            % Fuzzy Self-Tuning PSO
            d = delta(dim, pos(i, :), gBest);

            f = ( delta(dim, pos(i, :), pPosL(i, :)) / dmax ) * ( (min(fitness(i), fdelta) - min(pScoreL(i), fdelta) ) / abs(fdelta) );
            
            if ( 0 <= d < 0.2*dmax )
                Same = 1;
                Near = 0;
                Far = 0;
            elseif ( 0.2*dmax <= d < 0.4*dmax )
                Same = (0.4*dmax - d)/(0.2*dmax);
                Near = (d - 0.2*dmax)/(0.2*dmax);
                Far = 0;
            elseif ( 0.4*dmax <= d < 0.6*dmax )
                Same = 0;
                Near = (0.6*dmax - d)/(0.2*dmax);
                Far = (d - 0.4*dmax)/(0.2*dmax);
            else
                Same = 0;
                Near = 0;
                Far = 1;
            end
            
            if (f == -1)
                Better = 1;
                Worse = 0;
            elseif (-1 < f < 0)
                Better = -f;
                Worse = 0;
            elseif ((0 <= f) && (f < 1))
                Better = 0;
                Worse = f;
            else
                Better = 0;
                Worse = 1;
            end
            Unvaried = 1 - abs(f);
            
            w = ( (Worse + Same) * 0.3 + (Unvaried + Near) * 0.5 + (Better + Far) * 1 ) / (Worse + Same + Unvaried + Near + Better + Far);
            c2 = ( (Better + Near) * 1 + (Unvaried + Same) * 2 + (Worse + Far) * 3 ) / (Worse + Same + Unvaried + Near + Better + Far);
            c1 = ( (Far) * 0.1 + (Worse + Unvaried + Same + Near) * 1.5 + (Better) * 3 ) / (Worse + Same + Unvaried + Near + Better + Far);
            L = ( (Unvaried + Better + Far) * 0 + (Near + Same) * 0.001 + (Worse) * 0.01 ) / (Worse + Same + Unvaried + Near + Better + Far);
            U = ( (Same) * 0.1 + (Near + Better + Unvaried) * 0.15 + (Worse + Far) * 0.2 ) / (Worse + Same + Unvaried + Near + Better + Far);
            
            Vmax = U*(ub - lb);
            Vmin = L*(ub - lb);
        end
        
        r1 = rand();
        r2 = rand();
        for j=1:size(pos,2)       
            vel(i,j)=w*vel(i,j)+c1*r1*(pBest(i,j)-pos(i,j))+c2*r2*(gBest(j)-pos(i,j));
            
            if(vel(i,j)>Vmax)
                vel(i,j)=Vmax;
            end
            if(vel(i,j)<-Vmax)
                vel(i,j)=-Vmax;
            end
            pPosL(i, j) = pos(i, j);
            pos(i,j)=pos(i,j)+vel(i,j);
        end
        % Return back the particles that go beyond the boundaries of the search
        % space
        Flag4ub=pos(i,:)>ub;
        Flag4lb=pos(i,:)<lb;
        pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
    end
    
    cg_curve(l)=gBestScore;
	l=l+1;
    bestPos=gBest;
end

end

function [d] = delta(dim, pos1, pos2)
    tmp = 0;
    for j = 1 : dim
        tmp = tmp + ( pos1(j) - pos2(j) ) ^ 2;
    end
    d = sqrt(tmp);
end