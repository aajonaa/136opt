%_________________________________________________________________________%
%  Fick's Law Algorithm (FLA) source codes version 1.0                    %
%                                                                         %
%  Developed in MATLAB R2021b                                             %
%                                                                         %
%  Coresponding Author:  Abdelazim G. Hussien                             %
%                                                                         %
%                                                                         %
%         e-Mail: abdelazim.hussien@liu.se                                %
%                 aga08@fayoum.edu.eg                                     %
%                                                                         %
%                                                                         %
%   Main paper: Fatma Hashim, Reham R Mostafa, Abdelazim G. Hussien,      %
%                     Seyedali Mirjalili, & Karam M. Sallam               %
%               Knowledge-based Systems                                   %
%                                                                         %
%_________________________________________________________________________%
function [best_pos,  Convergence_curve] = FLA(N, Max_FEs, lb,ub, dim,fobj)
FEs = 0;
C1=0.5;C2=2;c3=.1;c4=.2;c5=2;
D=.01;
X=lb+rand(N,dim)*(ub-lb);%intial postions
AllFitness = inf * ones(1, N);
for i=1:N
    AllFitness(i) = fobj(X(i,:));
    FEs =  FEs + 1;
end
[bestFitness, IndexBestF] = min(AllFitness);
best_pos = X(IndexBestF,:);
n1=round(N/2);
n2=N-n1;
Fitness1 = inf * ones(1, n1);
Fitness2 = inf * ones(1, n2);
X1=X(1:n1,:);
X2=X(n1+1:N,:);
for i=1:n1
    Fitness1(i) = fobj(X1(i,:));
    FEs = FEs + 1;
end
for i=1:n2
    Fitness2(i) = fobj(X2(i,:));
    FEs = FEs + 1;
end
[Fit1_best, IndexFSeo1] = min(Fitness1);
[Fit2_best, IndexFSeo2] = min(Fitness2);
X1_best = X1(IndexFSeo1,:);
X2_best = X2(IndexFSeo2,:);
vec_flag=[1,-1];
if Fit1_best<Fit2_best
    FSss=Fit1_best;
    YSol=X1_best;
else
    FSss=Fit2_best;
    YSol=X2_best;
end

it = 1;
Convergence_curve = [];

% for FEs = 1:Max_FEs
while FEs <= Max_FEs
    TF(FEs)=sinh(FEs/Max_FEs)^C1;
    X=[X1;X2];
    %             DO
    if TF(FEs)<0.9
           DOF=exp(-(C2*TF(FEs)-rand))^C2;
           %DOF = exp(-(2*sinh(FEs/Max_FEs)^0.5-rand))^2
            TDO=c5*TF(FEs)-rand;%   direction of flow
            if (TDO)<rand
                %         select no of molecules
                M1N = c3*n1;
                M2N = c4*n1;
                NT12 =round((M2N-M1N).*rand(1,1) + M1N);
                for u=1:NT12
                    flag_index = floor(2*rand()+1);
                    DFg=vec_flag(flag_index);
                    Xm2=mean(X2);
                    Xm1=mean(X1);
                    J=-D*(Xm2-Xm1)/norm(X2_best- X1(u,:)+eps);
                    X1new(u,:)= X2_best+ DFg*DOF.*rand(1,dim).*(J.*X2_best-X1(u,:));
                end
                for u=NT12+1:n1
                    for tt=1:dim
                        p=rand;
                        if p<0.8
                            X1new(u,tt) = X1_best(tt);
                        elseif p<.9
                            r3=rand;
                            X1new(u,tt)=X1(u,tt)+DOF.*((ub-lb)*r3+lb);
                        else
                            X1new(u,tt) =X1(u,tt);
                        end
                        
                    end
                end
                for u=1:n2
                    r4=rand;
                    X2new(u,:)= X2_best+DOF.*((ub-lb)*r4+lb);
                end
            else
                M1N = .1*n2;
                M2N = .2*n2;
                Ntransfer =round((M2N-M1N).*rand(1,1) + M1N);
                for u=1:Ntransfer
                    flag_index = floor(2*rand()+1);
                    DFg=vec_flag(flag_index);
                    R1=randi(n1);
                    Xm1=mean(X1);
                    Xm2=mean(X2);
                    J=-D*(Xm1-Xm2)/norm(X1_best- X2(u,:)+eps);
                    X2new(u,:)=  X1_best+DFg*DOF.*rand(1,dim).*(J.*X1_best-1*X2(u,:));
                end
                for u=Ntransfer+1:n2
                    for tt=1:dim
                        p=rand;
                        if p<0.8
                            X2new(u,tt) = X2_best(tt);
                        elseif p<.9
                            r3=rand;
                            X2new(u,tt)=X2(u,tt)+DOF.*((ub-lb)*r3+lb);
                        else
                            X2new(u,tt) =X2(u,tt);
                        end
                        
                    end
                end
                for u=1:n1
                    r4=rand;
                    X1new(u,:)= X1_best+DOF.*((ub-lb)*r4+lb);
                end
            end
 
    else
%         Equilibrium operator (EO)
        if TF(FEs)<=1
            for u=1:n1
                flag_index = floor(2*rand()+1);
                DFg=vec_flag(flag_index);
                Xm1=mean(X1);
                Xmeo1=X1_best;
                J=-D*(Xmeo1-Xm1)/norm(X1_best- X1(u,:)+eps);
                DRF= exp(-J/TF(FEs));
                MS=exp(-Fit1_best/(Fitness1(u)+eps));
                R1=rand(1,dim);
                Qeo=DFg*DRF.*R1;
                X1new(u,:)= X1_best+Qeo.*X1(u,:)+Qeo.*(MS*X1_best-X1(u,:));
            end
            for u=1:n2
                flag_index = floor(2*rand()+1);
                DFg=vec_flag(flag_index);
                Xm2=mean(X2);
                Xmeo2=X2_best;
                J=-D*(Xmeo2-Xm2)/norm(X2_best- X2(u,:)+eps);
                DRF= exp(-J/TF(FEs));
                MS=exp(-Fit2_best/(Fitness2(u)+eps));
                R1=rand(1,dim);
                Qeo=DFg*DRF.*R1;
                X2new(u,:)=  X2_best+Qeo.*X2(u,:)+Qeo.*(MS*X1_best-X2(u,:));
            end
        else
            %     Steady state operator (SSO):
        for u=1:n1
            flag_index = floor(2*rand()+1);
            DFg=vec_flag(flag_index);
            Xm1=mean(X1);
            Xm=mean(X);
            J=-D*(Xm-Xm1)/norm(best_pos- X1(u,:)+eps);
            DRF= exp(-J/TF(FEs));
            MS=exp(-FSss/(Fitness1(u)+eps));
            R1=rand(1,dim);
            Qg=DFg*DRF.*R1;
            X1new(u,:)=  best_pos+Qg.*X1(u,:)+Qg.*(MS*best_pos-X1(u,:));
        end
        for u=1:n2
            Xm1=mean(X1);
            Xm=mean(X);
            J=-D*(Xm1-Xm)/norm(best_pos- X2(u,:)+eps);
            DRF= exp(-J/TF(FEs));
            MS=exp(-FSss/(Fitness2(u)+eps));
            flag_index = floor(2*rand()+1);
            DFg=vec_flag(flag_index);
                        Qg=DFg*DRF.*R1;
            X2new(u,:)= best_pos+ Qg.*X2(u,:)+Qg.*(MS*best_pos-X2(u,:));
        end
        end
    end
    for j=1:n1
        FU=X1new(j,:)>ub;FL=X1new(j,:)<lb;X1new(j,:)=(X1new(j,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        v = fobj(X1new(j,:));
        FEs = FEs + 1;
        if v<Fitness1(j)
            Fitness1(j)=v;
            X1(j,:)= X1new(j,:);
        end
    end
    for j=1:n2
        FU=X2new(j,:)>ub;FL=X2new(j,:)<lb;X2new(j,:)=(X2new(j,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        v = fobj(X2new(j,:));
        FEs = FEs + 1;
        if v<Fitness2(j)
            Fitness2(j)=v;
            X2(j,:)= X2new(j,:);
        end
    end
    
    [Fit1_best, IndexFSeo1] = min(Fitness1);
    [Fit2_best, IndexFSeo2] = min(Fitness2);
    
    X1_best = X1(IndexFSeo1,:);
    X2_best = X2(IndexFSeo2,:);
    if Fit1_best<Fit2_best
        FSss=Fit1_best;
        YSol=X1_best;
    else
        FSss=Fit2_best;
        YSol=X2_best;
    end
    Convergence_curve(it)=FSss;
    it = it + 1;
    if FSss<bestFitness
        bestFitness=FSss;
        best_pos =YSol;
        
    end 
end
end