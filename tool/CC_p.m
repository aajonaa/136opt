function [X,fitness,FEs] = CC_p(X,fitness,dim,lb,ub,fobj,p2)
    %% 
    % Criss
    FEs=0;
    Mhc = zeros(size(X,1),dim);
    Bhc = randperm(size(X,1));
    for i=1:(size(X,1)/2)
        no1= Bhc(2*i-1);
        no2 = Bhc(2*i);
        for j=1:dim
            r1 = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
            r2 = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
            c1 = (rand(1)*2)-1; %生成服从均匀分布的-1到1的随机数
            c2 = (rand(1)*2)-1;
            Mhc(no1,j)=r1*X(no1,j)+(1-r1)*X(no2,j)+c1*(X(no1,j)-X(no2,j));
            Mhc(no2,j)=r2*X(no2,j)+(1-r2)*X(no1,j)+c2*(X(no2,j)-X(no1,j));
        end
    end
     for i=1:size(X,1)
        % Check boundries
        FU=Mhc(i,:)>ub;FL=Mhc(i,:)<lb;Mhc(i,:)=(Mhc(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness_mhc(i)=fobj(Mhc(i,:));
        FEs=FEs+1;
        if fitness(i)<fitness_mhc(i)
            X(i,:)=X(i,:);
        else
            X(i,:)=Mhc(i,:);
            fitness(i)=fitness_mhc(i);
        end
     end    
     %Cross
    Bvc = randperm(dim);
    Mvc = X;
    %normalization
    for i=1:size(X,1)
        Boundary_no= size(ub,2); % numnber of boundaries
        % If the boundaries of all variables are equal and user enter a signle
        % number for both ub and lb
        if Boundary_no==1
            Mvc(i,:) = (Mvc(i,:)-lb)/(ub-lb);
        end
        % If each variable has a different lb and ub
        if Boundary_no>1
            for j=1:dim
                ub_j=ub(j);
                lb_j=lb(j);
                Mvc(i,j)=(Mvc(i,j)-lb_j)/(ub_j-lb_j);
            end
        end
    end
%     p2 = 0.6;  %p2取0.2到0.8之间
    for i=1:(dim/2)
        p = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
        if  p<p2
            no1= Bvc(2*i-1);
            no2 = Bvc(2*i);
            for j=1:size(X,1)
                r = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
                Mvc(j,no1)=r*Mvc(j,no1)+(1-r)*Mvc(j,no2);
            end
        end
    end
    
    for i=1:size(X,1)
        Boundary_no= size(ub,2); % numnber of boundaries
        % If the boundaries of all variables are equal and user enter a signle
        % number for both ub and lb
        if Boundary_no==1
            Mvc(i,:) = Mvc(i,:)*(ub-lb)+lb;
        end
        % If each variable has a different lb and ub
        if Boundary_no>1
            for j=1:dim
                ub_j=ub(j);
                lb_j=lb(j);
                Mvc(i,j)=(ub_j-lb_j)*Mvc(i,j)+lb_j;
            end
        end
        % Check boundries
        FU=Mvc(i,:)>ub;FL=Mvc(i,:)<lb;Mvc(i,:)=(Mvc(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness_mvc(i)=fobj(Mvc(i,:));
        FEs=FEs+1;
        if fitness(i)<fitness_mvc(i)
            X(i,:)=X(i,:);
        else
            X(i,:)=Mvc(i,:);
            fitness(i)=fitness_mvc(i);
        end
     end
    %%
end

