%**************************************************************************************************
%  Chen, Wei-Neng, et al. "Particle swarm optimization with an aging leader and challengers."
%  IEEE Transactions on Evolutionary Computation 17.2 (2013): 241-258.
%  ALC-PSO
%  Writer: Chen Xu
%  Date: 2017/8/9
%**************************************************************************************************


function [Leader,convergence]=ALCPSO(SearchAgents_no,maxFES,lb,ub,dim,fobj)
tic
if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
lu=[lb;ub];

rand('seed', sum(100 * clock));


D = dim;
Xmin = lu(1,:);
Xmax = lu(2,:);

popsize = SearchAgents_no;
w = 0.4;
c1 = 2; c2 = 2;
Age = 0;
lifespan = 60;
T = 2;
pro = 1/D;

% Step 1 Initialize the main population
X = repmat(Xmin, popsize, 1) + rand(popsize, D) .* (repmat(Xmax-Xmin, popsize, 1));
val_X=zeros(size(X,1),1);
for i=1:size(X,1)
    val_X(i) = fobj(X(i,:));
end
pBest = X; val_pBest = val_X;
[~,indexG] = min(val_pBest);
gBest = pBest(indexG,:); val_gBest = val_pBest(indexG,:);
Leader = gBest; val_Leader = val_gBest;

Vmax = (Xmax - Xmin)*0.5;    % velocity limit
Vmin = -Vmax;
V = zeros(popsize,D);

gen = 0;
FES = 0; 
convergence = [];  % record the best results
outcome_lifespan = lifespan*ones(popsize,1);
l=1;
while FES < maxFES
   
    % == == == == == Step 2 Velocity and Position Updating == == == == == %
    V = w*V + c1*rand(popsize,D).*(pBest-X) + c2*rand(popsize,D).*(repmat(Leader,popsize,1)-X);
    V = boundary_repair(V,Vmin,Vmax,'absorb');
    X = X + V ;
    X = boundary_repair(X,Xmin,Xmax,'reflect');
    
    % == == == == == Step 3 Updating pBest and Leader  == == == == == %
    val_X=zeros(size(X,1),1);
    for i=1:size(X,1)
        val_X(i) = fobj(X(i,:));
    end
    FES = FES+popsize;
    % update pBest
    previous_val_pBest = val_pBest;
    Index = (val_X<val_pBest);
    pBest(Index,:) = X(Index,:);  val_pBest(Index,:) = val_X(Index,:);
    % update Leader
    previous_val_Leader = val_Leader;
    [~,indexG] = min(val_X);
    if min(val_X)<val_Leader
        Leader = X(indexG(1),:);  val_Leader = val_X(indexG(1),:);
    end
    % update gBest
    previous_val_gBest = val_gBest;
    [~,indexG] = min(val_pBest);
    if min(val_pBest)<val_gBest
        gBest = pBest(indexG(1),:);  val_gBest = val_pBest(indexG(1),:);
    end
    
    % == == == == == Step 4 Lifespan Control   == == == == == %
    delta_val_pBest = sum(val_pBest - previous_val_pBest);
    delta_val_Leader = val_Leader - previous_val_Leader;
    delta_val_gBest = val_gBest - previous_val_gBest;
    if delta_val_gBest<0
        lifespan = lifespan + 2;
    elseif (delta_val_gBest==0) && (delta_val_pBest<0)
        lifespan = lifespan + 1;
    elseif (delta_val_pBest==0) && (delta_val_Leader<0)
        lifespan = lifespan;
    elseif delta_val_Leader==0
        lifespan = max(1,lifespan - 1);
    end
    
    Age = Age + 1;
    
    if Age>=lifespan
        % == == == == == Step 5 Generating a Challenger   == == == == == %
        Index = rand(1,D)<pro;
        if sum(Index)==0
            Index(randi(D)) = true;
        end
        Challenger = Leader;
        rr = rand(1,D);
        Challenger(Index) = Xmin(Index)+rr(Index).*(Xmax(Index)-Xmin(Index));
        val_Challenger = fobj(Challenger);
        FES = FES+1;
        if val_Challenger<val_gBest
            gBest = Challenger; val_gBest = val_Challenger;
        end
        
        % == == == == == Step 6 Evaluating the Challenger   == == == == == %
        delta_Ch = 0;
        for k = 1:T
            V = w*V + c1*rand(popsize,D).*(pBest-X) + c2*rand(popsize,D).*(repmat(Challenger,popsize,1)-X);
            V = boundary_repair(V,Vmin,Vmax,'absorb');
            X = X + V ;
            X = boundary_repair(X,Xmin,Xmax,'reflect');
            val_X=zeros(size(X,1),1);
            for i=1:size(X,1)
                val_X(i) = fobj(X(i,:));
            end
            FES = FES+popsize;
            % update pBest
            previous_val_pBest = val_pBest;
            Index = (val_X<val_pBest);
            pBest(Index,:) = X(Index,:);  val_pBest(Index,:) = val_X(Index,:);
            % update Challenger
            [~,indexG] = min(val_X);
            if min(val_X)<val_Challenger
                Challenger = X(indexG(1),:);  val_Challenger = val_X(indexG(1),:);
            end
            % update gBest
            [~,indexG] = min(val_pBest);
            if min(val_pBest)<val_gBest
                gBest = pBest(indexG(1),:);  val_gBest = val_pBest(indexG(1),:);
            end
            delta_Ch = delta_Ch+sum(previous_val_pBest-val_pBest);
        end
        
        if delta_Ch>0 % some pBest have been improved
            accept = 1;
            Leader = Challenger; Age = 0;
        else
            accept = 0;
            Age = Age-1;
        end
        
    end
%     convergence = [convergence val_gBest];
      convergence(l)=val_gBest;
      l=l+1;
    
end
%     bestScore=convergence(end);
toc
end






