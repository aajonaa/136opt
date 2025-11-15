function [ BestPosition,BestFitness,L,FES] = ADN_FES(Cr,Padn,L0,L,Lmin,Rou,dim,fobj)
%UNTITLED2 此处显示有关此函数的摘要
%   Cr待优化的最优解，Padn概率选择控制参数，Padn=0.2（原文默认值）,L0用于控制产生随机步长，L0=0.2（原文默认值）
%   L步长参数，Lmin最小步长，Rou收缩系数，dim问题维度，fobj适应度函数。
FES=0;
Nb=zeros(2*dim,dim);
Nb_fit=zeros(1,2*dim);
BestPosition=Cr;
BestFitness=fobj(Cr);
FES=FES+1;
if rand()<Padn
    for k=1:dim
        for d=1:dim
            if k==d
                Nb(2*k-1,d)=Cr(d)+L;
                Nb(2*k,d)=Cr(d)-L;
            else
                Nb(2*k-1,d)=Cr(d);
                Nb(2*k,d)=Cr(d);                
            end
        end
        Nb_fit(2*k-1)=fobj(Nb(2*k-1,:));
        Nb_fit(2*k)=fobj(Nb(2*k,:));
        FES=FES+2;
    end
    [bst_val,bst_id]=min(Nb_fit);
    if bst_val<BestFitness
        BestFitness=bst_val;
        BestPosition=Nb(bst_id,:);
    else
        L=L*Rou;%%收缩
        if L<Lmin
            L=L0*rand();%%扩张
        end
    end
end

end

