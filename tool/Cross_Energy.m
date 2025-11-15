%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%×ÝÏò½»²æ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitness_mvc,fitness_x,Mvc,X,Leader_pos,Leader_score,gbestfit,FEs,count]= Cross_Energy(X,Leader_pos,Leader_score,lb,ub,dim,FEs,fobj)
count = 0;
FU = X > ub;
FL = X < lb;
X = X .* (~(FU + FL)) + ub .* FU + lb .* FL;
Bvc = randperm(dim);
Mvc = X;
Boundary_size = size(ub,2);
if Boundary_size == 1
    Mvc = (Mvc - lb) / (ub - lb);
end
if Boundary_size > 1
    for j = 1:dim
        ub_j = ub(j);
        lb_j = lb(j);
        Mvc(j) = (Mvc(j) - lb_j) / (ub_j - lb_j);
    end
end
p2 = 1.0;
for j = 1:(dim/2)
    p = unifrnd(0,1);
    if p < p2
        no1 = Bvc(2 * j - 1);
        no2 = Bvc(2 * j);
        r = unifrnd(0,1);
        Mvc(no1) = r * Mvc(no1) + (1 - r) * Mvc(no2);
    end
end
if Boundary_size == 1
    Mvc = Mvc * (ub -lb) + lb;
end
if Boundary_size > 1
    for j = 1:dim
        ub_j = ub(j);
        lb_j = lb(j);
        Mvc(j) = (ub_j - lb_j) * Mvc(j) + lb_j;
    end    
end
FU = Mvc > ub;
FL = Mvc < lb;
Mvc = Mvc .* (~(FU + FL)) + ub .* FU + lb .* FL;
FU = X > ub;
FL = X < lb;
X = X .* (~(FU + FL)) + ub .* FU + lb .* FL;
fitness_x = fobj(X);
FEs = FEs + 1;
if fitness_x < Leader_score
    Leader_score = fitness_x;
    Leader_pos = X;
end
count = count + 1;
gbestfit(count) = Leader_score;
fitness_mvc = fobj(Mvc);
FEs = FEs + 1;
if fitness_mvc < Leader_score
    Leader_score = fitness_mvc;
    Leader_pos = Mvc;
end
count = count + 1;
gbestfit(count) = Leader_score;
end