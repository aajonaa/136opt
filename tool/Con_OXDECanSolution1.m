% function [ BestPosition,BestFitness] = Con_OXDECanSolution( Xi,gBest,dim,Q,F,fobj )
%% =================================this is for energy====================================
function [ BestPosition,BestFitness,fes] = Con_OXDECanSolution1( Xi,gBest,gbestfitness,dim,Q,F,fobj )

%%Xi个体i，gBest最优个体，dim问题维度，Q正交表水平数，F因素数，fobj适应度函数
%%BestPosition正交学习产生的最优解，BestFitness正交学习产生的最优解适应度函数值。
[ L,M ] = Get_QDOA( Q,F );
BestPosition=gBest;
% BestFitness=inf;
BestFitness=gbestfitness;
gbestscore=zeros(1,M+1);
fes=0;
count=0;
UB=max(Xi(1,:),gBest(1,:));
LB=min(Xi(1,:),gBest(1,:));
%% --------------------问题维度大于因素数，需要对问题维度分解分组映射到因素数上。
if dim>F
    
%% 产生F-1个随机整数序列，然后排序，用于问题维度分解映射到因素数上。
    a=2:dim-1;
    b=randperm(length(a));
    K=a(b(1:F-1));
    K=sort(K);
	%% 因素水平分解映射。
    for i=1:M
        for j=1:F
            if j==1
                Position(i,1:K(j))=LB(1:K(j))+(L(i,j)-1)*(UB(1:K(j))-LB(1:K(j)))./(Q-1);
            elseif j==F
                Position(i,K(j-1)+1:dim)=LB(K(j-1)+1:dim)+(L(i,j)-1)*(UB(K(j-1)+1:dim)-LB(K(j-1)+1:dim))./(Q-1);
            else
                Position(i,K(j-1)+1:K(j))=LB(K(j-1)+1:K(j))+(L(i,j)-1)*(UB(K(j-1)+1:K(j))-LB(K(j-1)+1:K(j)))./(Q-1);
            end
        end
        Fitness(1,i)=fobj( Position(i,:));
        fes=fes+1;
    end
	%% 因素水平分析构造最优水平因素组合
    for q=1:Q
        for i=1:M
            for j=1:F
                if L(i,j)==q
                    Z(i,j)=1;
                else
                    Z(i,j)=0;
                end
            end
        end
        S(q,:)=(Fitness*Z)./sum(Z); % S是一个QxF的矩阵 表示各个因素的个水平的影响因子
    end
    [~,Index]=max(S);
    for j=1:F
        if j==1
            Position(M+1,1:K(j))=LB(1:K(j))+(Index(j)-1)*(UB(1:K(j))-LB(1:K(j)))./(Q-1);
        elseif j==F
            Position(M+1,K(j-1)+1:dim)=LB(K(j-1)+1:dim)+(Index(j)-1)*(UB(K(j-1)+1:dim)-LB(K(j-1)+1:dim))./(Q-1);
        else
            Position(M+1,K(j-1)+1:K(j))=LB(K(j-1)+1:K(j))+(Index(j)-1)*(UB(K(j-1)+1:K(j))-LB(K(j-1)+1:K(j)))./(Q-1);
        end
    end
    Fitness(1,M+1)=fobj(Position(M+1,:));
    fes=fes+1;
	%% 计算所有候选解的适应度函数值，进而选出最优解
    for i=1:M+1
        if Fitness(1,i)<BestFitness
            BestFitness=Fitness(1,i);
            BestPosition=Position(i,:);
        end
        count=count+1;
        gbestscore(count)=BestFitness;
    end
    
%% -------------问题维度小于等于因素数，可以直接使用正交表
else 
    for i=1:M
        for j=1:dim
            Position(i,j)=LB(j)+(L(i,j)-1)*(UB(j)-LB(j))./(Q-1);
        end
        Fitness(1,i)=fobj( Position(i,:));
        fes=fes+1;
    end
    for q=1:Q
        for i=1:M
            for j=1:dim
                if L(i,j)==q
                    Z(i,j)=1;
                else
                    Z(i,j)=0;
                end
            end
        end
        S(q,:)=(Fitness*Z)./sum(Z);
    end
    [~,Index]=max(S);
    for j=1:dim
        Position(M+1,j)=LB(j)+(Index(j)-1)*(UB(j)-LB(j))./(Q-1);
    end
    Fitness(1,M+1)=fobj(Position(M+1,:));
    fes=fes+1;
    for i=1:M+1
        if Fitness(1,i)<BestFitness
            BestFitness=Fitness(1,i);
            BestPosition=Position(i,:);
        end
        count=count+1;
        gbestscore(count)=BestFitness;
    end
end
% FES_NEW=FES;
end
