function [best_fitness,chrom]=testQgate_mfo_test(chrom,fitness,best_fitness,best_pos,fobj)
%% 量子旋转门
sizepop=size(chrom,1);
lenchrom=size(chrom,2);
S=chrom;
for i=1:sizepop %Search_no
    for j=1:lenchrom %dim
        A=chrom(i,j);   % α
        B=best_pos(1,j);     % β
        delta=0;
        s=0;
        if (fitness(i)==best_fitness)
            delta=0;                  % delta为旋转角的大小
            s=0;                        % s为旋转角的符号，即旋转方向
        elseif (fitness(i)<best_fitness)
            delta=0.015*pi;
            if A*B<=0
                s=1;
            elseif A*B>0
                s=-1;
            elseif A==0
                s=0;
            elseif B==0
                s=sign(randn);
            end
        elseif (fitness(i)>best_fitness)
            delta=0.015*pi;
            if A*B<=0
                s=-1;
            elseif A*B>0
                s=1;
            elseif A==0
                s=sign(randn);
            elseif B==0
                s=0;
            end
        end
        e=s*delta;       % e为旋转角
        U=[cos(e) -sin(e);sin(e) cos(e)];      % 量子旋转门
        y=U*[A B]';        % y为更新后的量子位
        chrom(i,j)=y(1);    
    end
    Fitness1=fobj(chrom(i,:));
    Fitness2=fobj(S(i,:));
    Food_fitness=[Fitness1,Fitness2];
    [fit,m]=min(Food_fitness);
        if m==1   
            chrom(i,:)=chrom(i,:);
        else 
            chrom(i,:)=S(i,:); 
        end
        if fit<best_fitness
            best_fitness=fit;
        end
end
