%**************************************************************************************************
%  GOTLBO: Generalized Oppositional teaching learning based optimization
%  Writer: Chen Xu
%  Date: 2015/09/24
%**************************************************************************************************


function [GBEST, cg_curve]=GOTLBO(N,MaxFEs,lb,ub,dim,fobj)

if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
lu=[lb;ub];

rand('seed', sum(100 * clock));

GBEST = zeros(1,dim);
fes = 0;
D = dim;
n = D;
Xmin = lu(1,:);
Xmax = lu(2,:);

popsize = N;
Jr = 0.3;
%  Generalized Oppositional Population Initialization
X = repmat(Xmin, popsize, 1) + rand(popsize, n) .* (repmat(Xmax-Xmin, popsize, 1));
k = rand;
GOX = k*repmat(Xmin+Xmax, popsize, 1) - X;
GOX = boundary_repair(GOX, Xmin, Xmax,'random');
X = [X; GOX];
val_X=zeros(size(X,1),1);
for i=1:size(X,1)
  val_X(i) = fobj(X(i,:));
  fes = fes + 1;
end

% Select the NP fittest individuals
[~,index] = sort(val_X);
X = X(index,:); val_X = val_X(index,:);
X = X(1:popsize,:); val_X = val_X(1:popsize,:);


% FES = 2*popsize; 
% maxFES = D * 10000;
% convergence = [];  % record the best results
cg_curve=[];
l=1;
% while   FES < maxFES
while fes<MaxFEs    
    % == == == == == = Teacher Phase == == == == == =
    for i=1:popsize
        
        % Calculate the mean of each design variable
        mean_result = mean(X);
        % Identify the best solution (teacher)
        [val_Best,index_Best] = min(val_X);
        Best = X(index_Best(1),:);
        % Modify solution based on best solution
        TF=round(1+rand*(1));
        Xi = X(i,:) + (Best -TF*mean_result).*rand(1,D);
        Xi = boundary_repair(Xi,Xmin,Xmax,'reflect');
        % Accept or Reject
        val_Xi = fobj(Xi);
        fes = fes + 1;
%         FES = FES + 1;
        if val_Xi<val_X(i,:)
            val_X(i,:) = val_Xi; X(i,:) = Xi;
        end
        
    end
    
    % == == == == == = Student Phase == == == == == =
    for i=1:popsize
        
        j = randi(popsize);
        while j == i
            j = randi(popsize);
        end
        if val_X(i,:)<val_X(j,:)
            Xi = X(i,:) + rand(1,D).*(X(i,:)-X(j,:));
        else
            Xi = X(i,:) + rand(1,D).*(X(j,:)-X(i,:));
        end
        Xi = boundary_repair(Xi,Xmin,Xmax,'reflect');
        %  Accept or Reject
        val_Xi = fobj(Xi);
        fes = fes + 1;
%         FES = FES + 1;
        if  val_Xi<val_X(i,:)
            val_X(i,:) = val_Xi; X(i,:) = Xi;
        end
        
    end
    
    %  == == == == == = Generalized Opposition-Based Generation Jumping  == == == == == =
    if rand<Jr
        k = rand;
        GOX = k*repmat(max(X)+min(X),popsize,1) - X;
        GOX = boundary_repair_GOBL(GOX, Xmin, Xmax,min(X),max(X));
        val_GOX=zeros(size(GOX,1),1);
        for k=1:size(GOX,1)
          val_GOX(k) = fobj(GOX(k,:));
          fes = fes + 1;
        end
        
%         FES = FES + popsize;
        X = [X;GOX]; val_X = [val_X;val_GOX];
        % Select the NP fittest individuals
        [~,index] = sort(val_X);
        X = X(index,:); val_X = val_X(index,:);
        X = X(1:popsize,:); val_X = val_X(1:popsize,:);
    end
    
%     convergence = [convergence min(val_X)];
      cg_curve(l)=min(val_X);
      l=l+1;        
end

% bestScore=convergence(end);

end


