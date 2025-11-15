%**************************************************************************************************
%  BLPSO:Biogeography-Based Learning Particle Swarm Optimization
%  Writer: Xu Chen
%**************************************************************************************************

function [convergence]=BLPSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj,funcNum)
fobj = str2func(['cec14_func_', num2str(funcNum)]);

if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
lu=[lb;ub];

rand('seed', sum(100 * clock));


D = dim;


Xmin = lu(1,:);
Xmax = lu(2,:);

% BLPSO parameters
popsize = SearchAgents_no ;
maxFES = 1e4*D; maxGEN = maxFES/popsize;
iwt = 0.9 - (1 : maxGEN) * (0.7 / maxGEN);
c = 1.49445;
% BBO parameters
I = 1; % max immigration rate
E = 1; % max emigration rate

% Compute migration rates, assuming the population is sorted from most fit to least fit
MigrateModel = 5;
migration_models;

% Initialize the main population
X = repmat(Xmin, popsize, 1) + rand(popsize,D) .* (repmat(Xmax-Xmin, popsize, 1));
val_X=zeros(size(X,1),1);
for i=1:size(X,1)
    val_X(i) = fobj(X(i,:));
end
pBest = X; val_pBest = val_X;
[~,indexG] = min(val_pBest);
gBest = pBest(indexG,:); val_gBest = val_pBest(indexG,:);
Vmax = (Xmax - Xmin)*0.2;  Vmin = -Vmax;
V = repmat(Vmin,popsize,1) + rand(popsize,D).*repmat(Vmax-Vmin,popsize,1);


% FES = 0; 
GEN = 1;
% convergence = [];  % record the best results
convergence = zeros(1,Max_iteration);
l=1;
% while   FES < maxFES
while l< Max_iteration+1  
    for i = 1:popsize
        
        %  Biogeography-based exemplar generation method
        pBest_ind(i,:) = LearnIndex_BLPSO(val_pBest,popsize,D,i,mumu,lambda);
        
        %  Biogeography-based Learning Strategy
        for  j=1:D
            pBest_f(i,j) = pBest(pBest_ind(i,j),j);
        end
        V(i,:) = iwt(GEN)*V(i,:) + c*rand(1,D).*(pBest_f(i,:)-X(i,:));  % update velocity
        V(i,:) = boundConstraint_absorb(V(i,:),Vmin,Vmax);
        X(i,:) = X(i,:)+V(i,:);    % update position
        
        if all(X(i,:)<=Xmax) && all(X(i,:)>=Xmin)  % X(i,:) is feasible
            val = fobj(X(i,:));
%             FES = FES+1;
            if val<val_pBest(i)    % update pBest
                pBest(i,:) = X(i,:);  val_pBest(i) = val;
                if  val<val_gBest  % update gBest
                    gBest = X(i,:);  val_gBest = val;
                end
            end
        end
        
    end
    
    
    
%     convergence = [convergence val_gBest];
      convergence(l)=val_gBest;
      l=l+1;    
    GEN = GEN+1;
    if (GEN == maxGEN) && (l < Max_iteration+1)
        GEN = GEN-1;
    end
    
%     bestScore=convergence(end);
    
end





