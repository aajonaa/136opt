function [Best_rime,Convergence_curve]=LRRIME(N,MaxFEs,lb,ub,dim,fobj)
% RIME enhanced with learning + reflection strategy inspired by GO algorithm.
% title = {Growth Optimizer: A powerful metaheuristic algorithm for solving continuous and discrete global optimization problems},
% journal = {Knowledge-Based Systems},
% doi = {https://doi.org/10.1016/j.knosys.2022.110206},
% author = {Qingke Zhang and Hao Gao and Zhi-Hui Zhan and Junqing Li and Huaxiang Zhang},

Best_rime=zeros(1,dim);
Best_rime_rate=inf;
Rimepop=initialization(N,dim,ub,lb);
Lb=lb.*ones(1,dim); Ub=ub.*ones(1,dim);
FEs=0; Time = 1; Convergence_curve=[];
Rime_rates=zeros(1,N); newRime_rates=zeros(1,N);
W = 5; noImp = 0;

% eval init
for i=1:N
    Rime_rates(1,i)=fobj(Rimepop(i,:)); FEs=FEs+1;
    if Rime_rates(1,i)<Best_rime_rate
        Best_rime_rate=Rime_rates(1,i); Best_rime=Rimepop(i,:);
    end
end

while FEs < MaxFEs
    RimeFactor = (rand-0.5)*2*cos((pi*FEs/(MaxFEs/10)))*(1-round(FEs*W/MaxFEs)/W);
    E = sqrt(FEs/MaxFEs);
    newRimepop = Rimepop;
    normalized_rime_rates = normr(Rime_rates);

    % Soft-rime only (puncture disabled)
    for i=1:N
        for j=1:dim
            if rand()<E
                newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));
            end
        end
    end

    % Accept RIME moves greedily
    imp=false;
    for i=1:N
        Flag4ub=newRimepop(i,:)>ub; Flag4lb=newRimepop(i,:)<lb;
        xi=(newRimepop(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        nf=fobj(xi); FEs=FEs+1;
        if nf<Rime_rates(1,i)
            Rime_rates(1,i)=nf; Rimepop(i,:)=xi; imp=true;
            if nf< Best_rime_rate, Best_rime_rate=nf; Best_rime=xi; end
        end
        if FEs>=MaxFEs, Convergence_curve(Time)=Best_rime_rate; return; end
    end
    noImp = imp*0 + (~imp)*(noImp+1);

    % ---- GO miniblocks (learn + reflect) on elites ----
    [~, ind]=sort(Rime_rates,'ascend');
    ke = max(2, round(0.10*N)); elite = ind(1:ke);
    topK = min(5,N);
    mut = 0.01 + (0.10-0.01)*(1 - FEs/MaxFEs);

    % Option: gate on stagnation
    % if noImp >= 20
    % GO learning on half elites
    L = elite(1:max(1,round(0.4*numel(elite))));
    Best = Rimepop(ind(1),:);
    for ii=1:numel(L)
        i=L(ii);
        Worst=Rimepop(ind(randi([max(1,N-4),N])),:);
        Better=Rimepop(ind(randi([2,min(5,N)])),:);
        % two random distinct others (safe for small N)
        v = [1:i-1, i+1:N]; if numel(v) < 2, continue; end
        r1=v(randi(numel(v))); v(v==r1)=[]; r2=v(randi(numel(v)));
        D1=Best-Better; D2=Best-Worst; D3=Better-Worst; D4=Rimepop(r1,:)-Rimepop(r2,:);
        s = norm(D1)+norm(D2)+norm(D3)+norm(D4); if s==0, s=eps; end
        Gap = (norm(D1)/s)*D1 + (norm(D2)/s)*D2 + (norm(D3)/s)*D3 + (norm(D4)/s)*D4;
        SF = (Rime_rates(i)/max(Rime_rates));
        cand = Rimepop(i,:) + SF * Gap;
        cand = max(min(cand,ub),lb);
        if FEs>=MaxFEs, break; end
        nf=fobj(cand); FEs=FEs+1;
        if nf<Rime_rates(1,i)
            Rime_rates(1,i)=nf; Rimepop(i,:)=cand; imp=true;
            if nf< Best_rime_rate, Best_rime_rate=nf; Best_rime=cand; end
        end
        if FEs>=MaxFEs, break; end
    end

    % GO reflection on remaining elites
    R = elite(max(1, numel(elite)-max(1,round(0.4*numel(elite)))+1):end);
    for ii=1:numel(R)
        i=R(ii); cand=Rimepop(i,:);
        for j=1:dim
            if rand<0.30
                Rv=Rimepop(ind(randi(topK)),:);
                cand(j)=Rimepop(i,j) + (Rv(j)-Rimepop(i,j))*rand;
                if rand<mut, cand(j)=lb + (ub-lb).*rand(1,1); end
            end
        end
        cand=max(min(cand,ub),lb);
        if FEs>=MaxFEs, break; end
        nf=fobj(cand); FEs=FEs+1;
        if nf<Rime_rates(1,i)
            Rime_rates(1,i)=nf; Rimepop(i,:)=cand; imp=true;
            if nf< Best_rime_rate, Best_rime_rate=nf; Best_rime=cand; end
        end
        if FEs>=MaxFEs, break; end
    end
    % end % stagnation gating

    Convergence_curve(Time)=Best_rime_rate; Time=Time+1;
end
end

