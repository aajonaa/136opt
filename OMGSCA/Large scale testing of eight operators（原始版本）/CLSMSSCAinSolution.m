% 总述：
%
% CLSMSSCA：
% SCA + 多种群 MS（DNS、SRS、PDS） + 混沌局部搜索 CLS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 详述：
%
% 多种群（MS）：
% 由三个策略（DNS、SRS、PDS）组成。
% 论文引自：A multi-swarm particle swarm optimization algorithm based on dynamical topology and purposeful detecting。
% 
% 混沌局部搜索（CLS）：
% 
function [cg_curve]=CLSMSSCA(N,Max_iteration,lb,ub,dim,fobj,funcNum)
fobj = str2func(['cec14_func_', num2str(funcNum)]);
% Parameters of PSO
Vmax=6;
Vmin=-6;
noP=N;
wMax=0.9;
wMin=0.2;
c1=2;
c2=2;

% Parameters of MSPSO
iter=Max_iteration;
MaxFEs = dim * 10000;
fes = 0;
t = 0;
Rn = 10;
Stag = 0;
%NsubVector = [20, 10, 8, 5, 4, 2, 1];
j = 1;
for i = 2 : N / 2
    if ( mod(N/i, 1) == 0 )
        NsubVector(j) = N / i;
        j = j + 1;
    end
end
NsubVector(j) = 1;
m = 1;
Cgen = floor( iter/size(NsubVector, 2) );

% Divide each dimension search space into Rn same sized sub-regions
interval = ( ub - lb ) / Rn;
NumOfSuperior = round(Rn/3);
NumOfInferior = NumOfSuperior;
NumOfModerate = Rn - NumOfSuperior - NumOfInferior;

Nsub = NsubVector(m);
Ssub = N/Nsub;
M = zeros(Rn, dim);
tabu = zeros(Rn, dim);

% Initializations
pos=initialization(noP,dim,ub,lb);
SubSwarmIndex = randperm(N);
Swarm = zeros(Nsub, Ssub, dim);
k = 1;
for i = 1 : Nsub
    for j = 1 : Ssub
        Swarm(i, j, :) = pos(SubSwarmIndex(k), :);
        k = k + 1;
    end
end

vel=rand(Nsub, Ssub, dim) * (Vmax - Vmin) + Vmin;
pBestScore=zeros(Nsub, Ssub);
for i = 1 : Nsub
    for j = 1 : Ssub
        pBestScore(i, j)=inf;
    end
end
pBest=zeros(Nsub, Ssub, dim);
gBestScore = zeros(Nsub, 1);
for i = 1 : Nsub
    gBestScore(i) = inf;
end
gBest = zeros(Nsub, dim);
GBESTSCORE = inf;
cg_curve=zeros(1,iter);

% First time calculating the fitness of all particles
for i = 1 : Nsub
    for j = 1 : Ssub
        x = zeros(1, dim);
        for p = 1 : dim
            x(p) = Swarm(i, j, p);
        end
        fitness=fobj(x);

        if(pBestScore(i, j) > fitness)
            pBestScore(i, j) = fitness;
            pBest(i, j, :) = Swarm(i, j, :);
            if(gBestScore(i) > fitness)
                gBestScore(i) = fitness;
                gBest(i, :) = Swarm(i, j, :);
            end
        end
    end
end
[GBESTSCORE, gBestIndex] = min(gBestScore);
GBEST = gBest(gBestIndex,:);
minGBestScoreL = GBESTSCORE;

    
while ( t < iter )
    t = t + 1;
    
    'CLSMSSCA'
    t
    
    % Update the W of PSO
    %w=wMax-l*((wMax-wMin)/iter);
    %w=1;
    
    % Eq. (3.4)
    a = 2;
    Max_iteration = Max_iteration;
    r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
    
    for i = 1 : Nsub    
        % Update the Velocity and Position of particles
        for j = 1 : Ssub
            for k = 1 : dim
                
                % Update r2, r3, and r4 for Eq. (3.3)
                r2=(2*pi)*rand();
                r3=2*rand;
                r4=rand();
                
                % Eq. (3.3)
                if r4<0.5
                    % Eq. (3.1)
                    Swarm(i,j,k)= Swarm(i,j,k)+(r1*sin(r2)*abs(r3*gBest(i,k)-Swarm(i,j,k)));
                else
                    % Eq. (3.2)
                    Swarm(i,j,k)= Swarm(i,j,k)+(r1*cos(r2)*abs(r3*gBest(i,k)-Swarm(i,j,k)));
                end
                
            end
            %% Chaos
                K = N;
                % start chaos
                setCan = (Max_iteration-t+1)/Max_iteration;
                x = rand();
                while(~(x~= 0.25 && x~= 0.5 && x~=0.75))
                    x = rand();
                end
                ch(1) = x;
                for ant=1:(K)
                    ch(ant+1)=4*ch(ant)*(1-ch(ant));
                    CH(ant,:) = lb+ch(ant)*(ub-lb);    %ub大
                    V = (1-setCan)*gBest(i, :)+setCan*CH(ant);
                   % ObjValV=feval(objfun,V);         %计算函数值
                    Flag4ub=V>ub;
                    Flag4lb=V<lb;
                    V=(V.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                    FitnessV=fobj(V);%计算适应度值
                    %    if (ObjValV<GlobalMin)
                    if (FitnessV<gBestScore(i))
                        gBestScore(i) = FitnessV;
                        gBest(i, :) = V;
                        %break;
                    end
                end
                % end chaos
                %%
        end
        
        for j = 1 : Ssub
            % Return back the particles that go beyond the boundaries of the search
            % space
            Flag4ub=Swarm(i, j, :)>ub;
            Flag4lb=Swarm(i, j, :)<lb;
            Swarm(i, j, :)=(Swarm(i, j, :).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            
            % Calculate objective function for each particle
            x = zeros(1, dim);
            for p = 1 : dim
                x(p) = Swarm(i, j, p);
            end
            fitness=fobj(x);
            fes = fes + 1;

            if(pBestScore(i, j) > fitness)
                pBestScore(i, j) = fitness;
                pBest(i, j, :) = Swarm(i, j, :);
                if(gBestScore(i) > fitness)
                    gBestScore(i) = fitness;
                    gBest(i, :) = Swarm(i, j, :);
                end
            end
        end
    end
    
    % calculate the MER
    for i = 1 : Nsub
        for j = 1 : Ssub
            for k = 1 : dim
                p = ceil( (pBest(i, j, k) - lb) / interval );
                if ( p == 0 )
                    p = 1;
                end
                M(p, k) = M(p, k) + 1;
            end
        end
    end
    
    % update the GBEST and Stagbest
    if ( min(gBestScore) < GBESTSCORE )
        [GBESTSCORE, gBestIndex] = min(gBestScore);
        GBEST = gBest(gBestIndex,:);
    end
    if ( min(gBestScore) == minGBestScoreL )
        Stag = Stag + 1;
    else
        minGBestScoreL = min(gBestScore);
        Stag = 0;
    end
    
    % DNS
    [fes, Nsub, Ssub, m, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, gBestIndex] = DNS(fes, t, Cgen, N, Nsub, Ssub, m, NsubVector, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, fobj, dim, gBestIndex);
    
    % SRS
    [Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest] = SRS(Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, N, dim);
    
    % PDS
    [fes, GBEST, tabu, GBESTSCORE] = PDS(fes, GBEST, M, tabu, GBESTSCORE, NumOfSuperior, NumOfModerate, NumOfInferior, interval, dim, lb, fobj);
    
    cg_curve(t) = GBESTSCORE;
end

end

function [fes, Nsub, Ssub, m, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, gBestIndex] = DNS(fes, t, Cgen, N, Nsub, Ssub, m, NsubVector, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, fobj, dim, gBestIndex)
if ( mod(t, Cgen) == 0 && m < size(NsubVector, 2) )
    m = m + 1;
    NsubL = Nsub;
    Nsub = NsubVector(m);
    SsubL = Ssub;
    Ssub = N/Nsub;
    
    % Random divide the whole population into Nsub sub-swarms
    SubSwarmIndex = randperm(N);
    SwarmL = Swarm;
    Swarm = zeros(Nsub, Ssub, dim);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), SsubL);
            if (tmp == 0)
                tmp = SsubL;
            end
            Swarm(i, j, :) = SwarmL(ceil(SubSwarmIndex(k)/SsubL), tmp, :);
            k = k + 1;
        end
    end
    
    velL = vel;
    vel = zeros(Nsub, Ssub, dim);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), SsubL);
            if (tmp == 0)
                tmp = SsubL;
            end
            vel(i, j, :) = velL(ceil(SubSwarmIndex(k)/SsubL), tmp, :);
            k = k + 1;
        end
    end
    
    pBestScoreL = pBestScore;
    pBestScore = zeros(Nsub, Ssub);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), SsubL);
            if (tmp == 0)
                tmp = SsubL;
            end
            pBestScore(i, j) = pBestScoreL(ceil(SubSwarmIndex(k)/SsubL), tmp);
            k = k + 1;
        end
    end
    
    pBestL = pBest;
    pBest = zeros(Nsub, Ssub, dim);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), SsubL);
            if (tmp == 0)
                tmp = SsubL;
            end
            pBest(i, j, :) = pBestL(ceil(SubSwarmIndex(k)/SsubL), tmp, :);
            k = k + 1;
        end
    end
    
    gBestScore = zeros(Nsub, 1);
    TmpIndex = zeros(Nsub, 1);
    gBest = zeros(Nsub, dim);
    for i = 1 : Nsub
        [gBestScore(i), TmpIndex(i)] = min(pBestScore(i, :));
        gBest(i, :) = pBest(i, TmpIndex(i), :);
    end
    
    % to carry out BFGS Quasi-Newton method
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'MaxIter', floor(0.1*fes), 'MaxFunEvals', floor(0.1*fes));
    [x, fval]  = fminunc(fobj, GBEST, options);
    if (fval < GBESTSCORE)
        GBESTSCORE = fval;
        GBEST = x;
    end
    fes = fes + 0.1*fes;
end
end

function [Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest] = SRS(Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, N, dim)
if (Stag >= Ssub/2)
    SubSwarmIndex = randperm(N);
    SwarmL = Swarm;
    Swarm = zeros(Nsub, Ssub, dim);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), Ssub);
            if (tmp == 0)
                tmp = Ssub;
            end
            Swarm(i, j, :) = SwarmL(ceil(SubSwarmIndex(k)/Ssub), tmp, :);
            k = k + 1;
        end
    end
    
    velL = vel;
    vel = zeros(Nsub, Ssub, dim);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), Ssub);
            if (tmp == 0)
                tmp = Ssub;
            end
            vel(i, j, :) = velL(ceil(SubSwarmIndex(k)/Ssub), tmp, :);
            k = k + 1;
        end
    end
    
    pBestScoreL = pBestScore;
    pBestScore = zeros(Nsub, Ssub);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), Ssub);
            if (tmp == 0)
                tmp = Ssub;
            end
            pBestScore(i, j) = pBestScoreL(ceil(SubSwarmIndex(k)/Ssub), tmp);
            k = k + 1;
        end
    end
    
    pBestL = pBest;
    pBest = zeros(Nsub, Ssub, dim);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), Ssub);
            if (tmp == 0)
                tmp = Ssub;
            end
            pBest(i, j, :) = pBestL(ceil(SubSwarmIndex(k)/Ssub), tmp, :);
            k = k + 1;
        end
    end
    
    gBestScore = zeros(Nsub, 1);
    gBestIndex = zeros(Nsub, 1);
    gBest = zeros(Nsub, dim);
    for i = 1 : Nsub
        [gBestScore(i), gBestIndex(i)] = min(pBestScore(i, :));
        gBest(i, :) = pBest(i, gBestIndex(i), :);
    end
    
    Stag = 0;
end
end

function [fes, GBEST, tabu, GBESTSCORE] = PDS(fes, GBEST, M, tabu, GBESTSCORE, NumOfSuperior, NumOfModerate, NumOfInferior, interval, dim, lb, fobj)
[~, MorTabuIndex] = sort(M, 1, 'DESCEND');

TmpGB = GBEST;
for i = 1 : dim
    Tmp = [];
    k = ceil( (TmpGB(i) - lb) / interval );
    if (k == 0)
        k = 1;
    end
    p = find(MorTabuIndex(:, i) == k);
    if (p <= NumOfSuperior)
        IsInferiorTabuFull = 1;
        
        for j = NumOfSuperior + NumOfModerate + 1 : NumOfSuperior + NumOfModerate + NumOfInferior
            if (tabu(MorTabuIndex(j, i), i) == 0)
               Tmp = [Tmp; lb + interval*(MorTabuIndex(j, i)-1 + rand())];
               IsInferiorTabuFull = 0;
            end
        end
        
        if (IsInferiorTabuFull == 0)
            X = randperm(size(Tmp, 1));
            TmpGB(i) = Tmp(X(1));
            
            fes = fes + 1;
            TmpGBScore = fobj(TmpGB);
            if (TmpGBScore < GBESTSCORE)
              GBEST = TmpGB;
              GBESTSCORE = TmpGBScore;
            end

            tabu(MorTabuIndex(j, i), i) = 1;
        end
        
        if (any(tabu(:, i) == 0) == 0)
            tabu(:, i) = tabu(:, i) * 0;
        end
    end
end
end