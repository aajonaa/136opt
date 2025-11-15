% %%%%%%%%%%%%%%%%生成拉普拉斯的随机数实现拉普拉斯交叉算法%%%%%%%%%%%%%%%%%%%%
% mu=0;                      %均值
% 
% sigma=1;                  %标准差，方差的开平方
% 
% b=sigma/sqrt(2);      %根据标准差求相应的b
% 
% a=rand(1,10000)-0.5;    %生成(-0.5,0.5)区间内均匀分布的随机数列 (一万个数的行向量);
% 
% x=mu-b*sign(a).*log(1-2*abs(a)); %生成符合拉普拉斯分布的随机数列
function [Leader_pos,Leader_score,gbestfit,fes,count]= Lapulas1(X,Leader_pos,Leader_score,lb,ub,dim,fes,belt,fobj)
count=0;
alph=0;
Mhc1 = zeros(1,dim);
Mhc2 = zeros(1,dim);
% worst_score=inf;
% worst_pos=zeros(1,dim);

% [fitnessX,index] = sort(fitnessX,'ascend');
% worst_score = fitnessX(N);  
% worst_pos = X(N,:);
rIndexArray = randperm(size(X,1));

      for j =1:dim
%         alph=Laplacian(0);
%         belt=Laplacian(0);
        u=rand(1);
        if u<=0.5                             %%拉普拉斯公式
           k=alph-belt*log(u);                   
        else
           k=alph+belt*log(u);
        end       
        Mhc1(j) = X(rIndexArray(1),j) + k * abs(X(rIndexArray(1),j) - Leader_pos(j)); 
        Mhc2(j) = Leader_pos(j) + k * abs(X(rIndexArray(1),j) - Leader_pos(j));
      end


    FU =  Mhc1 > ub;
    FL =  Mhc1 < lb;
    Mhc1 = (Mhc1 .* (~(FU + FL))) + ub .* FU + lb .* FL;
    fitness_mhc1 = fobj(Mhc1);
    fes = fes + 1;
    count = count + 1;
    
    if fitness_mhc1 < Leader_score
        Leader_score = fitness_mhc1;
        Leader_pos = Mhc1;
    end
    
    gbestfit(count) = Leader_score;

    FU =  Mhc2 > ub;
    FL =  Mhc2 < lb;
    Mhc2 = (Mhc2 .* (~(FU + FL))) + ub .* FU + lb .* FL;
    fitness_mhc2 = fobj(Mhc2);
    fes = fes + 1;
    count = count + 1;

    if fitness_mhc2 < Leader_score
        Leader_score = fitness_mhc2;
        Leader_pos = Mhc2;
    end
    gbestfit(count) = Leader_score;

    
    
    
    
end


% %%%%%%%%%%%%%%%%%%%将拉普拉斯分布函数写成函数形式%%%%%%%%%%%%%%%%%%
% function L=Laplacian(mu)
% sigma = 1;
% b=sigma/sqrt(2);
% a=rand()-0.5;
% L=mu-b*sign(a).*log(1-2*abs(a));
% end
