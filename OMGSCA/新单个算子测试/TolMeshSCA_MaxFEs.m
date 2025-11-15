% 总述：
%
% OLMSSCA：
% SCA + 多种群 + 正交学习 OL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 详述：
%
% 多种群（MS）：
% 由三个策略（DNS、SRS、PDS）组成。
% 论文引自：A multi-swarm particle swarm optimization algorithm based on dynamical topology and purposeful detecting
% 
function [GBEST,cg_curve]=SCA_MaxFEs(N,MaxFEs,lb,ub,dim,fobj)
if (size(ub, 2)~=1)
    dim = size(ub, 2);
else
    lb = lb*ones(1,dim);
    ub = ub*ones(1,dim);
end
% Parameters of PSO
noP=N;

% Parameters of MSPSO
iter=1000;
fes = 0;
t = 0;
Stag = 0;

Nsub = 1;
Ssub = N;

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
cg_curve=[];

% Parameters for OLMSSCA2
Objective_values = zeros(Nsub, Ssub);

% First time calculating the fitness of all particles
for i = 1 : Nsub
    for j = 1 : Ssub
        x = zeros(1, dim);
        for p = 1 : dim
            x(p) = Swarm(i, j, p);
        end
        fitness=fobj(x);
        fes = fes + 1;
        Objective_values(i, j) = fitness;

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


    
while ( fes < MaxFEs )
    t = t + 1;
    
    'SCA_MaxFEs'
    t
    
    % Update the W of PSO
    %w=wMax-l*((wMax-wMin)/iter);
    %w=1;
    
    % Eq. (3.4)
    a = 2;
    r1=a-fes*((a)/MaxFEs); % r1 decreases linearly from a to 0
    
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
        end
        
        for j = 1 : Ssub
            % Return back the particles that go beyond the boundaries of the search
            % space
            x = zeros(1, dim);
            for k = 1 : dim
                Flag4ub(k)=Swarm(i, j, k)>ub(k);
                Flag4lb(k)=Swarm(i, j, k)<lb(k);
                x(k) = Swarm(i, j, k);
            end
            x=(x.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            Swarm(i, j, :) = x;
            
            % Calculate objective function for each particle
            x = zeros(1, dim);
            for p = 1 : dim
                x(p) = Swarm(i, j, p);
            end
            fitness=fobj(x);
            Objective_values(i, j) = fitness;
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
        
        % type3
        poptions=psoptimset('TolFun',1e-15,'TolX', 1e-15,'TolMesh',1e-15,'MaxFunEvals', floor(0.1*fes));
        [outX, outf, ~, output] = patternsearch(fobj,gBest(i, :),[],[],[],[],[],[],[],poptions);

        if ( outf < gBestScore(i) )
            gBestScore(i) = outf;
            gBest(i, :) = outX;
        end
        fes = fes + output.funccount;
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
    %[fes, Nsub, Ssub, m, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, gBestIndex, Trial, Objective_values, V] = DNS(fes, t, Cgen, N, Nsub, Ssub, m, NsubVector, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, fobj, dim, gBestIndex, Trial, Objective_values, V);
    %VFitness=zeros(Nsub, Ssub);
    
    % SRS
    %[Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, Trial, Objective_values, V] = SRS(Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, N, dim, Trial, Objective_values, V);
    
    % PDS
    %[fes, GBEST, tabu, GBESTSCORE] = PDS(fes, GBEST, M, tabu, GBESTSCORE, NumOfSuperior, NumOfModerate, NumOfInferior, interval, dim, lb, fobj);
    
    cg_curve(t) = GBESTSCORE;
end

end

function [fes, Nsub, Ssub, m, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, gBestIndex, Trial, Objective_values, V] = DNS(fes, t, Cgen, N, Nsub, Ssub, m, NsubVector, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, fobj, dim, gBestIndex, Trial, Objective_values, V)
if ( mod(t, Cgen) == 0 && m < size(NsubVector, 2) )
    m = m + 1;
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
    
    TrialL = Trial;
    Trial = zeros(Nsub, Ssub);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), SsubL);
            if (tmp == 0)
                tmp = SsubL;
            end
            Trial(i, j) = TrialL(ceil(SubSwarmIndex(k)/SsubL), tmp);
            k = k + 1;
        end
    end
    
    Objective_valuesL = Objective_values;
    Objective_values = zeros(Nsub, Ssub);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), SsubL);
            if (tmp == 0)
                tmp = SsubL;
            end
            Objective_values(i, j) = Objective_valuesL(ceil(SubSwarmIndex(k)/SsubL), tmp);
            k = k + 1;
        end
    end
    
    VL = V;
    V = zeros(Nsub, Ssub, dim);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), SsubL);
            if (tmp == 0)
                tmp = SsubL;
            end
            V(i, j, :) = VL(ceil(SubSwarmIndex(k)/SsubL), tmp, :);
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

function [Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, Trial, Objective_values, V] = SRS(Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, N, dim, Trial, Objective_values, V)
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
    
    TrialL = Trial;
    Trial = zeros(Nsub, Ssub);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), Ssub);
            if (tmp == 0)
                tmp = Ssub;
            end
            Trial(i, j) = TrialL(ceil(SubSwarmIndex(k)/Ssub), tmp);
            k = k + 1;
        end
    end
    
    Objective_valuesL = Objective_values;
    Objective_values = zeros(Nsub, Ssub);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), Ssub);
            if (tmp == 0)
                tmp = Ssub;
            end
            Objective_values(i, j) = Objective_valuesL(ceil(SubSwarmIndex(k)/Ssub), tmp);
            k = k + 1;
        end
    end
    
    VL = V;
    V = zeros(Nsub, Ssub, dim);
    k = 1;
    for i = 1 : Nsub
        for j = 1 : Ssub
            tmp = mod(SubSwarmIndex(k), Ssub);
            if (tmp == 0)
                tmp = Ssub;
            end
            V(i, j, :) = VL(ceil(SubSwarmIndex(k)/Ssub), tmp, :);
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
    k = ceil( (TmpGB(i) - lb(i)) / interval(i) );
    if (k == 0)
        k = 1;
    end
    p = find(MorTabuIndex(:, i) == k);
    if (p <= NumOfSuperior)
        IsInferiorTabuFull = 1;
        
        for j = NumOfSuperior + NumOfModerate + 1 : NumOfSuperior + NumOfModerate + NumOfInferior
            if (tabu(MorTabuIndex(j, i), i) == 0)
               Tmp = [Tmp; lb(i) + interval*(MorTabuIndex(j, i)-1 + rand())];
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