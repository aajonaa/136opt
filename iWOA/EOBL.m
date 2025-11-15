function [ P,fes ] = EOBL( X,dim,fobj )
%EOBL 此处显示有关此函数的摘要
%   此处显示详细说明
fes=0;

ub_oppo=max(X);
lb_oppo=min(X);

for i=1:size(X,1)
    fitness(i,:)=fobj(X(i,:));
    fes=fes+1;
end

for i=1:size(X,1)
    GOX(i,:)=rand()*(ub_oppo+lb_oppo)-X(i,:);
    for j=1:dim
        if GOX(i,j)>ub_oppo(j) || GOX(i,j)<lb_oppo(j)
            GOX(i,j)=rand()*(ub_oppo(j)-lb_oppo(j))+lb_oppo(j);
        end
    end
    Gfitness(i,:)=fobj(GOX(i,:));
    fes=fes+1;
end

allX=[X;GOX];
allfitness=[fitness;Gfitness];
[~,index]=sort(allfitness);
P=allX(index,:);
P=P(1:size(X,1),:);

end

