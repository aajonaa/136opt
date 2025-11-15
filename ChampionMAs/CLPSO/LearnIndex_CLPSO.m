function  pBest_ind = LearnIndex_CLPSO(pfit,popsize,D,i,Pci);
%  exemplar index

for j = 1:D
    if rand<Pci %产生的随机数小于pci，则该粒子学习其它粒子
        ind1 = randi(popsize);%产生一个1到popsize的随机整数
        while ind1==i, ind1 = randi(popsize);  end  %若产生的随机数等于i，则从新产生一个随机数
        ind2 = randi(popsize); 
        while ind2==i||ind2==ind1, ind2 = randi(popsize);  end
        if pfit(ind1)<pfit(ind2)
            pBest_ind(j) = ind1;
        else
            pBest_ind(j) = ind2;
        end 
    else
        pBest_ind(j) = i;
    end 
end

if all(pBest_ind==i)    
    ind = randi(popsize);
    while ind==i, ind = randi(popsize);  end
    pBest_ind(randi(D)) = ind; 
end

