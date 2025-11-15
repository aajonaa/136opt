function [next_par par_best]=select_and_recombine(par,par_best,n)
    F=0.5;          %加速因子
    CR=0.8;         %变异率
    next_par=par;   %新种群
    for i=1:n
        
        while 1         %在原种群中任选三个互不相同的个体进行交叉变异
            r1=floor(1+20*rand());
            r2=floor(1+20*rand());
            r3=floor(1+20*rand());
            if i~=r1 && i~=r2 && i~=r3 &&...
               r1~=r2 && r1~=r3 && r2~=r3
                break;
            end
        end
                     
        if rand()<CR    %变异率，可以对每一个特征分别设置，我这里要变一起变了
            next_par(i).x=par(r1).x+F*(par(r2).x-par(r3).x);    %交叉变异准则
            next_par(i).y=par(r1).y+F*(par(r2).y-par(r3).y);              
        end
        
        %计算变异后个体的适应度
        next_par(i).fit=compute_fit(next_par(i));
        %如果新个体没有变异前个体适应度高，新个体还原为旧个体
        if par(i).fit>next_par(i).fit
            next_par(i)=par(i);
        end
        %如果变异后适应度高于种群最高适应个体，则更新种群适应度最高个体
        if next_par(i).fit>par_best.fit
            par_best=next_par(i);
        end
    end    
end

function re=compute_fit(par)
    x=par.x;
    y=par.y;
    sigma=50;
    if x<-100 || x>100 || y<-100 || y>100
        re=0;        %超出范围适应度为0
    else            %否则适应度按目标函数求解
        re=(1/(2*pi*sigma^2))*exp(-(x.^2+y.^2)/(2*sigma^2)); 
    end
end