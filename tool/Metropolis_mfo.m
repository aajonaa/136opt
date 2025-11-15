function [R2,L,Fes]=Metropolis_mfo(S2,fobj,SearchAgents_no,dim,ub,lb,Food_pos,Fes)
%%Levy+混沌
    %S1:新解,行向量
    %S2：原来得解，行向量
    %dim
    %N:最大迭代次数
    %L：对应的行向量
    %% 参数设置
    K=SearchAgents_no;
    MaxIt=30;     % Maximum Number of Iterations

    MaxSubIt=10;    % Maximum Number of Sub-iterations

    T0=0.1;       % Initial Temp.

    alpha=0.99;     % Temp. Reduction Rate
    
    %% 数值
Fes = Fes +1;

     R2=fobj(S2);
     L=S2;
     T=T0;
     for i=1:SearchAgents_no
         Fes = Fes+2;

        S1=S2+rand(1,dim).*Levy(dim);
        S3=chaos(ub,lb,i,SearchAgents_no,Food_pos);
%         S3=S1*(1+randn(1));
        R1=fobj(S1);
        R3=fobj(S3);
        a=[R1,R3];
        [fit,m]=min(a);
        if m==2
            S1=S3;
        end
        if fit<R2
            S2=S1;
            R2=fit;
            L=S2;
        elseif exp(-(fit-R2)/T)>=rand
            S2=S1;
            R2=fit; 
            L=S2;
        else
            L=S2;
        end
        T=T*alpha;
     end
        
end
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end

