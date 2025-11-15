function  pBest_ind = LearnIndex_BLPSO(pfit,popsize,D,i,mu,lambda)
% Biogeography-based exemplar generation method

[~,sort_index] = sort(pfit);
for k = 1:length(sort_index)
    index(sort_index(k)) = k;
end 

mu1 = mu(index);
lambda1 = lambda(index); 

for k = 1:D
    if rand < lambda1(i)  %  Should we immigrate?
        % Yes - Pick a solution from which to emigrate (roulette wheel selection)
        RandomNum = rand * sum(mu1);
        Select = mu1(1);  
        SelectIndex = 1;
        while (RandomNum > Select) && (SelectIndex < popsize)
            SelectIndex = SelectIndex + 1;
            Select = Select + mu1(SelectIndex);
        end
        pBest_ind(k) =  SelectIndex; 
    else
        pBest_ind(k) = i;
    end 
end

if all(pBest_ind == i)    
    ind = randi(popsize);
    while ind==i, ind = randi(popsize);  end
    pBest_ind(randi(D)) = ind; 
end