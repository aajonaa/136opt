%**************************************************************************************************
%**************************************************************************************************

function [best_score,best_pos,time,Convergence_curve]=TLBO(SearchAgents_no,MaxFES,lb,ub,dim,fobj)
t = cputime;
%Initialize some parameters
best_pos=zeros(1,dim);
best_score=1; %change this to -inf for maximization problems

Convergence_curve=[];
FES=0;
fitcount=0;
% l=0;


if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
lu=[lb;ub];

rand('seed', sum(100 * clock));


D = dim;
n = D;
Xmin = lu(1,:);
Xmax = lu(2,:);

popsize = SearchAgents_no;
%  Generalized Oppositional Population Initialization
X = initialization(popsize,dim,ub,lb);   %随机初始化种群   %repmat(A,m,n)将矩阵A赋值m*n块
val_X=zeros(size(X,1),1);     %存放所有个体的适应度值
for i=1:size(X,1)
  val_X(i) = fobj(X(i,:)); 
  FES=FES+1;
  if val_X(i)<best_score
     best_score=val_X(i);
	 best_pos=X(i,:);
  end
  fitcount=fitcount+1;
  Convergence_curve(1,fitcount)=best_score;  
end

while   FES < MaxFES
%while l<Max_iteration+1    
    % == == == == == = Teacher Phase == == == == == =
%     FES
%     'TLBO'
    for i=1:popsize
        
        % Calculate the mean of each design variable
        mean_result = mean(X);     %计算种群个体的均值
        % Identify the best solution (teacher)
        [best_score,index_Best] = min(val_X);
        best_pos = X(index_Best(1),:);
        % Modify solution based on best solution
        TF=round(1+rand*(1));  %round对小数进行四舍五入运算的，这里TF是教学因素，被启发式的设置为1或者2，所以这里就看rand产生的数值是否大于0.5，若大于0.5则TF=2，若小于0.5，则TF=1
        Xi = X(i,:) + (best_pos -TF*mean_result).*rand(1,D); %教学阶段
        Xi = boundary_repair(Xi,Xmin,Xmax,'reflect');
        % Accept or Reject
        val_Xi = fobj(Xi);
        FES = FES + 1; 
        if val_Xi<val_X(i,:)
            val_X(i,:) = val_Xi;
            X(i,:) = Xi;
        end
        
        if val_X(i,:)<best_score
            best_score=val_X(i,:);
            best_pos=X(i,:);
        end
        fitcount=fitcount+1;
        Convergence_curve(1,fitcount)=best_score;
    end
    
    % == == == == == = Student Phase == == == == == =
    for i=1:popsize       
        j = randi(popsize);
        while j == i  %论文中是说随机选择两个个体，这里是将遍历到的个体做为其中一个，再选择一个与之不相等的个体        
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
        FES = FES + 1;
        if  val_Xi<val_X(i,:)
            val_X(i,:) = val_Xi;
            X(i,:) = Xi;
        end
        if val_X(i,:)<best_score
            best_score=val_X(i,:);
            best_pos=X(i,:);
        end
        fitcount=fitcount+1;
        Convergence_curve(1,fitcount)=best_score;
        
    end
end
Convergence_curve=Convergence_curve(:,(fitcount-MaxFES+1):end);
time = cputime - t;
end          


function u = boundary_repair(v,low,up,str)

[NP, D] = size(v);   
u = v; 
llow = repmat(low,NP,1);
uup  = repmat(up,NP,1);

if strcmp(str,'absorb')
    index1 = (u>uup);
    u(index1) = uup(index1); 
    index2 = (u<llow);
    u(index2) = llow(index2); 
end
   

if strcmp(str,'random')
    index1 = (u>uup)|(u<llow);
    rr = rand(size(u));
    u(index1) = llow(index1) + rr(index1).*(uup(index1)-llow(index1));  
end


if strcmp(str,'reflect')
    index1 = (u>uup);
    u(index1) = max( 2*uup(index1)-u(index1), llow(index1) );
    index2 = (u<llow);
    u(index2) = min( 2*llow(index2)-u(index2), uup(index2) );
end
end

