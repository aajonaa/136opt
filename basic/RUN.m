%-----------------------------评估框架----------------------------------------------------------------------------------------------

% RUNge Kutta optimizer (RUN)
% RUN Beyond the Metaphor: An Efficient Optimization Algorithm Based on Runge Kutta Method
% Website of RUN:http://imanahmadianfar.com/codes/
% Website of RUN:http://www.aliasgharheidari.com/RUN.html

% Iman Ahmadianfar : Assistant Professor, Department of Civil Engineering, Behbahan Khatam Alanbia University of Technology, Behbahan, Iran
% Ali asghar Heidari,Amir H. Gandomi , Xuefeng  Chu, Huiling Chen  

%  Thanks to Professor Amir H Gandomi (citations above 19000)
%  Faculty of Engineering & Information Technology, University of Technology Sydney, NSW 2007, Australia

%  Last update: 04-22-2021

%  e-Mail: im.ahmadian@gmail.com,i.ahmadianfar@bkatu.ac.ir.
%  e-Mail: as_heidari@ut.ac.ir, aliasghar68@gmail.com,
%  e-Mail (Singapore): aliasgha@comp.nus.edu.sg, t0917038@u.nus.edu
%---------------------------------------------------------------------------------------------------------------------------
%  Co-author: Ali Asghar Heidari(as_heidari@ut.ac.ir),Amir H Gandomi(a.h.gandomi@gmail.com),Xuefeng Chu(xuefeng.chu@ndsu.edu),Huiling Chen(chenhuiling.jlu@gmail.com), 
%---------------------------------------------------------------------------------------------------------------------------

% After use, please refer to the main paper:
% Iman Ahmadianfar, Ali Asghar Heidari,Amir H Gandomi,Xuefeng Chu,Huiling Chen,  
% N Beyond the Metaphor: An Efficient Optimization Algorithm Based on Runge Kutta Method
% Expert Systems With Applications, https://..... (Q1, 5-Year Impact Factor: 5.448, H-INDEX: 184)
%---------------------------------------------------------------------------------------------------------------------------
% Please follow the paper for related updates in researchgate: https://www.researchgate.net/profile/Iman_Ahmadianfar
%  Researchgate: https://www.researchgate.net/profile/Ali_Asghar_Heidari.

%  Website of RUN:%  http://www.aliasgharheidari.com/RUN.html,http://www.imanahmadianfar.com.

% You can also use and compare with our other new optimization methods:
                                                                       %(GBO)-2020-http://www.imanahmadianfar.com/codes.
                                                                       %(HGS)-2021- http://www.aliasgharheidari.com/HHO.html
                                                                       %(SMA)-2020- http://www.aliasgharheidari.com/SMA.html
                                                                       %(HHO)-2019- http://www.aliasgharheidari.com/HHO.html                                                                       


function [best_pos,Convergence_curve]=RUN(N,Max_FEs,lb,ub,dim,fobj)

tic

Cost=zeros(N,1);                % Record the Fitness of all Solutions
X=initialization(N,dim,ub,lb);  % Initialize the set of random solutions   
Xnew2=zeros(1,dim);

Convergence_curve = [];

% for i=1:N
%     Cost(i) = fobj(X(i,:));      % Calculate the Value of Objective Function
% end

FEs=0;
for i=1:N
    Cost(i) = fobj(X(i,:));     % Calculate the Value of Objective Function
    FEs=FEs+1;
end


[Best_Cost,ind] = min(Cost);     % Determine the Best Solution
best_pos = X(ind,:);

Convergence_curve = [];


it=1;%Number of iterations
%% Main Loop of RUN 

while FEs < Max_FEs
%     it=it+1;
    f=20.*exp(-(12.*(FEs/Max_FEs))); % (Eq.17.6) 
    Xavg = mean(X);               % Determine the Average of Solutions
    SF=2.*(0.5-rand(1,N)).*f;    % Determine the Adaptive Factor (Eq.17.5)
    
    for i=1:N
            [~,ind_l] = min(Cost);
            lBest = X(ind_l,:);   
            
            [A,B,C]=RndX(N,i);   % Determine Three Random Indices of Solutions
            [~,ind1] = min(Cost([A B C]));
            
            % Determine Delta X (Eqs. 11.1 to 11.3)
            gama = rand.*(X(i,:)-rand(1,dim).*(ub-lb)).*exp(-4*FEs/Max_FEs);  
            Stp=rand(1,dim).*((best_pos-rand.*Xavg)+gama);
            DelX = 2*rand(1,dim).*(abs(Stp));
            
            % Determine Xb and Xw for using in Runge Kutta method
            if Cost(i)<Cost(ind1)                
                Xb = X(i,:);
                Xw = X(ind1,:);
            else
                Xb = X(ind1,:);
                Xw = X(i,:);
            end

            SM = RungeKutta(Xb,Xw,DelX);   % Search Mechanism (SM) of RUN based on Runge Kutta Method
                        
            L=rand(1,dim)<0.5;
            Xc = L.*X(i,:)+(1-L).*X(A,:);  % (Eq. 17.3)
            Xm = L.*best_pos+(1-L).*lBest;   % (Eq. 17.4)
              
            vec=[1,-1];
            flag = floor(2*rand(1,dim)+1);
            r=vec(flag);                   % An Interger number 
            
            g = 2*rand;
            mu = 0.5+.1*randn(1,dim);
            
            % Determine New Solution Based on Runge Kutta Method (Eq.18) 
            if rand<0.5
                Xnew = (Xc+r.*SF(i).*g.*Xc) + SF(i).*(SM) + mu.*(Xm-Xc);
            else
                Xnew = (Xm+r.*SF(i).*g.*Xm) + SF(i).*(SM)+ mu.*(X(A,:)-X(B,:));
            end  
            
        % Check if solutions go outside the search space and bring them back
        FU=Xnew>ub;FL=Xnew<lb;Xnew=(Xnew.*(~(FU+FL)))+ub.*FU+lb.*FL; 
        CostNew=fobj(Xnew);
        FEs=FEs+1;
        if CostNew<Cost(i)
            X(i,:)=Xnew;
            Cost(i)=CostNew;
        end
%% Enhanced solution quality (ESQ)  (Eq. 19)      
        if rand<0.5
            EXP=exp(-5*rand*FEs/Max_FEs);
            r = floor(Unifrnd(-1,2,1,1));

            u=2*rand(1,dim); 
            w=Unifrnd(0,2,1,dim).*EXP;               %(Eq.19-1)
            
            [A,B,C]=RndX(N,i);
            Xavg=(X(A,:)+X(B,:)+X(C,:))/3;           %(Eq.19-2)         
            
            beta=rand(1,dim);
            Xnew1 = beta.*(best_pos)+(1-beta).*(Xavg); %(Eq.19-3)
            
            for j=1:dim
                if w(j)<1 
                    Xnew2(j) = Xnew1(j)+r*w(j)*abs((Xnew1(j)-Xavg(j))+randn);
                else
                    Xnew2(j) = (Xnew1(j)-Xavg(j))+r*w(j)*abs((u(j).*Xnew1(j)-Xavg(j))+randn);
                end
            end
            
            FU=Xnew2>ub;FL=Xnew2<lb;Xnew2=(Xnew2.*(~(FU+FL)))+ub.*FU+lb.*FL;
            CostNew=fobj(Xnew2);            
            FEs=FEs+1;
            if CostNew<Cost(i)
                X(i,:)=Xnew2;
                Cost(i)=CostNew;
            else
                if rand<w(randi(dim)) 
                    SM = RungeKutta(X(i,:),Xnew2,DelX);
                    Xnew = (Xnew2-rand.*Xnew2)+ SF(i)*(SM+(2*rand(1,dim).*best_pos-Xnew2));  % (Eq. 20)
                    
                    FU=Xnew>ub;FL=Xnew<lb;Xnew=(Xnew.*(~(FU+FL)))+ub.*FU+lb.*FL;
                    CostNew=fobj(Xnew);
                    FEs=FEs+1;
                    if CostNew<Cost(i)
                        X(i,:)=Xnew;
                        Cost(i)=CostNew;
                    end
                end
            end
        end
% End of ESQ         
%% Determine the Best Solution
        if Cost(i)<Best_Cost
            best_pos=X(i,:);
            Best_Cost=Cost(i);
        end

    end
% Save Best Solution at each iteration    
Convergence_curve(it) = Best_Cost;
it=it+1;
% disp(['it : ' num2str(it) ', Best Cost = ' num2str(Convergence_curve(it) )]);
% disp(['RUN iteration : ' num2str(it) ', Best Cost = ' num2str(Best_Cost )]);
end
toc
end

% A funtion to determine a random number 
%with uniform distribution (unifrnd function in Matlab) 
function z=Unifrnd(a,b,c,dim)
a2 = a/2;
b2 = b/2;
mu = a2+b2;
sig = b2-a2;
z = mu + sig .* (2*rand(c,dim)-1);
end

% A function to determine thress random indices of solutions
function [A,B,C]=RndX(nP,i)
Qi=randperm(nP);Qi(Qi==i)=[];
A=Qi(1);B=Qi(2);C=Qi(3);
end










