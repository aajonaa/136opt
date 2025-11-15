% 总述：
%
% COVMSSCA：
% SCA + 多种群 + 协方差 COV
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 详述：
%
% 多种群（MS）：
% 由三个策略（DNS、SRS、PDS）组成。
% 论文引自：A multi-swarm particle swarm optimization algorithm based on dynamical topology and purposeful detecting
% 
function [cg_curve]=COVMSSCA(N,Max_iteration,lb,ub,dim,fobj,funcNum)
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
%% 针对COVMSSCA的特殊处理
%因为Covariance这个函数需要size(RouBestPosition,1)大于等于2，即需要RouNum大于等于2，又因为RouNum =
%ceil(Ssub/4)，所以需要种群数目大于4
i = 0;
while ( N / NsubVector(i+1) <= 4 )
    i = i + 1;
end
if ( i ~= 0 )
    for j = 1 : i
        NsubVector(1) = [];
    end
end
%NsubVector(1) = [];
%%
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

%% 为协方差SCA做准备
% Parameters of COV
VFitness = zeros(Nsub, Ssub);
Objective_values = zeros(Nsub, Ssub);
VPosition = zeros(Nsub, Ssub, dim);
RouNum = ceil(Ssub/4);
RouBestPosition = zeros(RouNum, dim);
%%


% First time calculating the fitness of all particles
for i = 1 : Nsub
    for j = 1 : Ssub
        x = zeros(1, dim);
        for p = 1 : dim
            x(p) = Swarm(i, j, p);
        end
        fitness=fobj(x);
        %% 为协方差SCA做准备
        Objective_values(i, j) = fitness;
        %%

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
    
    'COVMSSCA'
    t
    
    % Update the W of PSO
    %w=wMax-l*((wMax-wMin)/iter);
    %w=1;
    
    % Eq. (3.4)
    a = 2;
    r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
    
    for i = 1 : Nsub
        %% 为协方差SCA做准备
        [~,Sort_Index]=sort(Objective_values(i, :));
        for q = 1 : RouNum
            RouBestPosition(q, :) = Swarm(i, Sort_Index(q), :);
        end
        %%
        
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
                    VPosition(i,j,k)= Swarm(i,j,k)+(r1*sin(r2)*abs(r3*gBest(i,k)-Swarm(i,j,k)));
                else
                    % Eq. (3.2)
                    VPosition(i,j,k)= Swarm(i,j,k)+(r1*cos(r2)*abs(r3*gBest(i,k)-Swarm(i,j,k)));
                end
                
                %% 为协方差SCA做准备
                x = zeros(1, dim);
                for p = 1 : dim
                    x(p) = VPosition(i, j, p);
                end
                VFitness(i, j) = fobj(x);
                
                if VFitness(i, j) < Objective_values(i, j)
                    Objective_values(i, j) = VFitness(i, j);
                    Swarm(i, j, :) = VPosition(i, j, :);
                end
            end
        end
        
        %% 协方差SCA
        for j = 1 : Ssub
            [ ZPosition ] = Covariance( RouBestPosition );
            if fobj(ZPosition) < Objective_values(i, j)
                Objective_values(i, j) = fobj(ZPosition);
                Swarm(i, j, :) = ZPosition;
            end
        end
        %% 协方差SCA结束
        
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
            %% 为协方差SCA做准备
            Objective_values(i, j) = fitness;
            %%
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
    
    %% 为协方差SCA做准备
    % DNS 会重塑整个多种群的划分，即会改变Nsub、Ssub这两个值，并随机重置
    [fes, Nsub, Ssub, m, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, gBestIndex, Objective_values] = DNS(fes, t, Cgen, N, Nsub, Ssub, m, NsubVector, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, fobj, dim, gBestIndex, Objective_values);
    
    % SRS 会随机重置整个多种群
    [Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, Objective_values] = SRS(Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, N, dim, Objective_values);
    VFitness = zeros(Nsub, Ssub);
    VPosition = zeros(Nsub, Ssub, dim);
    RouNum = ceil(Ssub/4);
    RouBestPosition = zeros(RouNum, dim);
    %%
    
    % PDS 仅仅针对GBESTSCORE进行优化
    [fes, GBEST, tabu, GBESTSCORE] = PDS(fes, GBEST, M, tabu, GBESTSCORE, NumOfSuperior, NumOfModerate, NumOfInferior, interval, dim, lb, fobj);
    
    cg_curve(t) = GBESTSCORE;
end

end

function [fes, Nsub, Ssub, m, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, gBestIndex, Objective_values] = DNS(fes, t, Cgen, N, Nsub, Ssub, m, NsubVector, Swarm, vel, pBestScore, pBest, gBestScore, gBest, GBESTSCORE, GBEST, fobj, dim, gBestIndex, Objective_values)
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
    
    %% 为协方差SCA做准备
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
    %%
    
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

function [Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, Objective_values] = SRS(Stag, Nsub, Ssub, Swarm, vel, pBestScore, pBest, gBestScore, gBest, N, dim, Objective_values)
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
    
    %% 为协方差SCA做准备
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
    %%
    
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