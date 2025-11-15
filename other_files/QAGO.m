%% Mod Jona 2025-1-8
function [best_pos,Convergence_curve]= QAGO(N,MaxFEs,lb,ub,dim,fobj)

%---------------------------------------------------------------------------
% Algorithm Name: QAGO
% gbestx: The global best solution ( gbestx = [x1,x2,...,xD]).
% gbestfitness: Record the fitness value of global best individual.
% gbesthistory: Record the history of changes in the global optimal fitness.
%---------------------------------------------------------------------------

%% (1)Initialization
FEs=0;
X=unifrnd(lb,ub,N,dim);
bestFitness=inf;
AllFitness=inf(N,1);
Convergence_curve = [];
it = 1;
for i=1:N
    AllFitness(i)=fobj(X(i,:));
    FEs=FEs+1;
    if bestFitness>=AllFitness(i)
        bestFitness=AllFitness(i);
        best_pos=X(i,:);
    end
    Convergence_curve(it)=bestFitness;
end

%% (2) Loop iteration
while FEs<=MaxFEs
    [~, ind]=sort(AllFitness);
    % Parameter adaptation based on distribution
    P1=ceil(unifrnd(0.05,0.2)*N);
    P2=normrnd(0.001*ones(1,N),0.001*ones(1,N));
    P3=normrnd(0.3*rand(N,dim),0.01);
    
    %% 1. Improved learning phase
    % Sampling individuals
    Best_X=X(ind(1),:);
    worse_index=ind(randi([N-P1+1,N],N,1));
    Worst_X=X(worse_index,:);
    better_index=ind(randi([2,P1],N,1));
    Better_X=X(better_index,:);
    normal_index=ind(randi([P1+1,N-P1],N,1));
    Normal_X=X(normal_index,:);
    [L1,L2,L3,L4]=selectID(N);
    
    for i=1:N
        Gap(1,:)=(Best_X-Better_X(i,:));
        Gap(2,:)=(Better_X(i,:)-Normal_X(i,:));
        Gap(3,:)=(Normal_X(i,:)-Worst_X(i,:));
        Gap(4,:)=(X(L1(i),:)-X(L2(i),:));
        Gap(5,:)=(X(L3(i),:)-X(L4(i),:));
        
        % Parameter self-adaptation based on one-dimensional mapping of vectors
        DGap(1,:)=(Best_X*Better_X(i,:)');
        DGap(2,:)=(Better_X(i,:)*Normal_X(i,:)');
        DGap(3,:)=(Normal_X(i,:)*Worst_X(i,:)');
        DGap(4,:)=(X(L1(i),:)*X(L2(i),:)');
        DGap(5,:)=(X(L3(i),:)*X(L4(i),:)');
        minDistance=min(DGap);
        DGap=DGap+2*abs(minDistance);
        LF=DGap./sum(DGap);
        
        % Parameter self-adaptation based on fitness difference
        FGap(1,:)=(abs(AllFitness(ind(1))-AllFitness(better_index(i))));
        FGap(2,:)=(abs(AllFitness(better_index(i))-AllFitness(normal_index(i))));
        FGap(3,:)=(abs(normal_index(i)-AllFitness(worse_index(i))));
        FGap(4,:)=(abs(AllFitness(L1(i))-AllFitness(L2(i))));
        FGap(5,:)=(abs(AllFitness(L3(i))-AllFitness(L4(i))));
        SF=FGap./sum(FGap);
        
        % Parameter self-adaptation based on Jensen-Shannon divergence
        LS=(LF+SF)/2;
        Djs=0.5*sum(LF.*log(LF./LS))+0.5*sum(SF.*log(SF./LS));
        djs=sqrt(Djs);
        
        % Learning operator refinement
        newX(i,:)=X(i,:)+sum(Gap.*(djs.*LF+(1-djs).*SF),1);
        
        % Boundary constraints
        Flag4ub= newX(i,:)>ub;
        Flag4lb= newX(i,:)<lb;
        newX(i,:)=(newX(i,:).*(~(Flag4ub+Flag4lb)))+(lb+(ub-lb)*rand(1,dim)).*Flag4ub+(lb+(ub-lb)*rand(1,dim)).*Flag4lb;
        
        % Evaluation
        newAllFitness(i)=fobj(newX(i,:));
        FEs=FEs+1;
        
        % Selection
        if AllFitness(i)>=newAllFitness(i)
            AllFitness(i)=newAllFitness(i);
            X(i,:)=newX(i,:);
            if bestFitness>=AllFitness(i)
                bestFitness=AllFitness(i);
                best_pos=X(i,:);
            end
        else
            if rand<P2(i)&&ind(i)~=ind(1)
                AllFitness(i)=newAllFitness(i);
                X(i,:)=newX(i,:);
            end
        end
    end % end for
    
    %% 2. Improved reflection phase
    newX=X;
    P2=normrnd(0.001*ones(1,N),0.001*ones(1,N));
    VSCR=rand(N,dim);
    VSAF=rand(N,dim);
    AF=0.01*(1-FEs/MaxFEs);
    I1=VSCR<P3;
    I2=VSAF<AF;
    R=ind(randi(P1,N,dim));
    
    for i=1:N
        % Reflection operator refinement
        for j=1:dim
            if I1(i,j)
                if I2(i,j)
                    newX(i,j)=lb+(ub-lb)*rand;
                else
                    S=randperm(N,3);
                    S(S==i)=[];
                    S(S==R(i,j))=[];
                    RM=S(1);
                    newX(i,j)=X(i,j)+rand*((X(R(i,j),j)-X(i,j))+(X(RM,j)-X(i,j)));
                end
            end
        end
        
        % Boundary constraints
        flagub=newX(i,:)>ub;
        newX(i,flagub)=(X(i,flagub)+ub)/2;
        flaglb=newX(i,:)<lb;
        newX(i,flaglb)=(X(i,flaglb)+lb)/2;
        
        % Evaluation
        newAllFitness(i)=fobj(newX(i,:));
        FEs=FEs+1;
        
        % Selection
        if AllFitness(i)>=newAllFitness(i)
            AllFitness(i)=newAllFitness(i);
            X(i,:)=newX(i,:);
            if bestFitness>=AllFitness(i)
                bestFitness=AllFitness(i);
                best_pos=X(i,:);
            end
        else
            if rand<P2(i)&&ind(i)~=ind(1)
                AllFitness(i)=newAllFitness(i);
                X(i,:)=newX(i,:);
            end
        end
    end % end for

    [~, idx] = sort(AllFitness);
    best_pos = X(idx, :);

    it = it + 1;
    Convergence_curve(it)=bestFitness;
    
    % fprintf("QAGO, FEs: %d, fitess error = %e\n",FEs,gbestfitness);
    if FEs>=MaxFEs
        break;
    end
    
end % end while

end % end QAGO algorithm
%--------------------------end-------------------------------------------------

%Subfunction:  Select four individuals different from Xi
function [L1,L2,L3,L4]=selectID(N)
for i=1:N
    if 2 <= N
        vecc=[1:i-1,i+1:N];
        r=zeros(1,4);
        for kkk =1:4
            n = N-kkk;
            t = randi(n,1,1);
            r(kkk) = vecc(t);
            vecc(t)=[];
        end
        L1(i)=r(1);L2(i)=r(2);L3(i)=r(3);L4(i)=r(4);
    end
end
L1=L1';L2=L2';L3=L3';L4=L4';
end