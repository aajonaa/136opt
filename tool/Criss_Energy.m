%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ºáÏò½»²æ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mhc,fitness_mhc,X,Leader_pos,Leader_score,gbestfit,FEs,count]= Criss_Energy(X,Leader_pos,Leader_score,lb,ub,dim,FEs,fobj)
count = 0;
Mhc = zeros(size(X,1),dim);
Bhc = randperm(size(X,1));
for i = 1:(size(X,1)/2)
    no1 = Bhc(2 * i -1);
    no2 = Bhc(2 * i);
    for j =1:dim 
        r1 = unifrnd(0,1);
        r2 = unifrnd(0,1);
        c1 = (rand(1) * 2) - 1;
        c2 = (rand(1) * 2) - 1;
        Mhc(no1,j) = r1 * X(no1,j) + (1 - r1) * X(no2,j) + c1 * (X(no1,j) - X(no2,j));
        Mhc(no2,j) = r2 * X(no2,j) + (1 - r2) * X(no1,j) + c2 * (X(no2,j) - X(no1,j));
    end
end
for i = 1:size(X,1)
    FU = Mhc(i,:) > ub;
    FL = Mhc(i,:) < lb;
    Mhc(i,:) = (Mhc(i,:) .* (~(FU + FL))) + ub .* FU + lb .* FL;
    fitness_mhc(i) = fobj(Mhc(i,:));
    FEs = FEs + 1;
    if fitness_mhc(i) < Leader_score
         Leader_score = fitness_mhc(i);
         Leader_pos = Mhc(i,:);
    end
    count = count + 1;
    gbestfit(count) = Leader_score;
end
    
end