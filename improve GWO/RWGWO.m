function [Alpha_pos,Convergence_curve]=RWGWO(SearchAgents_no,MaxFes,lb,ub,D,fobj)

Positions=initialization(SearchAgents_no,D,ub,lb);
Fes = 0;
for i=1:size(Positions,1)
    %fitness1(i)=fhd((Positions(i,:))',fobj);
    Fes = Fes + 1;
    fitness1(i) = fobj(Positions(i,:));
end

fitness=fitness1';

New_mat1=[fitness Positions];
[values, order]=sort(New_mat1(:,1));
New_mat2=New_mat1(order,:); % New_mat2第一列是排序后的适应度值 其余保存的是按照适应度值排序后的种群位置

Alpha_score= New_mat2(1,1); % Update alpha
Alpha_pos= New_mat2(1,2:size( New_mat2,2));  % Update alpha Score
Beta_pos= New_mat2(2,2:size( New_mat2,2));% Update Beta
Delta_pos= New_mat2(3,2:size( New_mat2,2));% Update Delta

Positions=New_mat2(:,2:size(New_mat2,2));
fitness=New_mat2(:,1);
U_Positions=zeros(SearchAgents_no,D);

it=1;
Convergence_curve=[];
while Fes<MaxFes
    xo=0;
    gamma=1;
    par=2-2*(Fes/MaxFes);
    for i=1:3
        for j=1:D
            y=rand();
            randomwalk=xo+gamma*tan(pi*(y-0.5));
            k=rand();
            U_Positions(i,j)=Positions(i,j)+randomwalk*par;
        end
    end
    
    
    a=2-2*(Fes/MaxFes);   % linearly decreasing vector from 2 to 0
    
    % Update the Position of search agents including omegas
    for i=4:size(Positions,1)
        for j=1:size(Positions,2)
            
            r1=rand();
            r2=rand();
            
            A1=2*a*r1-a;
            C1=2*r2;
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j));
            X1=Alpha_pos(j)-A1*D_alpha;
            
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a;
            C2=2*r2;
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j));
            X2=Beta_pos(j)-A2*D_beta;
            
            r1=rand();
            r2=rand();
            
            A3=2*a*r1-a;
            C3=2*r2;
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j));
            X3=Delta_pos(j)-A3*D_delta;
            
            U_Positions(i,j)=(X1+X2+X3)/3;
            
        end
        
    end
    
    for i=1:size(U_Positions,1)
%         for j=1:size(U_Positions,2)
%             if U_Positions(i,j)>ub
%                 U_Positions(i,j)=ub;
%             end
%             if U_Positions(i,j)<lb
%                 U_Positions(i,j)=lb;
%             end
%         end
        Flag4ub=U_Positions(i,:)>ub;
        Flag4lb=U_Positions(i,:)<lb;
        U_Positions(i,:)=(U_Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %U_fitness(i)=fhd((U_Positions(i,:))',fobj);
        Fes = Fes + 1;
        U_fitness(i)=fobj(U_Positions(i,:));
        if U_fitness(i)>fitness(i)
            U_fitness(i)=fitness(i);
            U_Positions(i,:)=Positions(i,:);
        end
    end
    
    
    
    New_mat1=[U_fitness' U_Positions];
    [values, order]=sort(New_mat1(:,1));
    New_mat2=New_mat1(order,:);
    Alpha_score= New_mat2(1,1); % Update alpha
    Alpha_pos= New_mat2(1,2:size( New_mat2,2));
    Beta_pos= New_mat2(2,2:size( New_mat2,2));
    Delta_pos= New_mat2(3,2:size( New_mat2,2));
    Positions=New_mat2(:,2:size(New_mat2,2));
    fitness=New_mat2(:,1);
    Convergence_curve(it)=Alpha_score;
    it=it+1;
end
end


