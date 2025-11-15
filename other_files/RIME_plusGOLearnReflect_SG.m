function [Best_rime,Convergence_curve]=RIME_plusGOLearnReflect_SG(N,MaxFEs,lb,ub,dim,fobj)
% RIME + GO learn+reflect, stagnation-gated, with light dedup pruning
% - Keep soft-rime; disable puncture
% - Trigger GO miniblocks only after no improvement for T iterations
% - Adaptive elite fraction and reflection prob with stagnation
% - No random acceptance; greedy only; mutation decays over FEs

Best_rime=zeros(1,dim);
Best_rime_rate=inf;
Rimepop=initialization(N,dim,ub,lb);
Lb=lb.*ones(1,dim); Ub=ub.*ones(1,dim);
FEs=0; Time = 1; Convergence_curve=[];
Rime_rates=zeros(1,N); newRime_rates=zeros(1,N);
W = 5; noImp = 0; lastBest = inf;

% init eval
for i=1:N
    Rime_rates(1,i)=fobj(Rimepop(i,:)); FEs=FEs+1;
    if Rime_rates(1,i)<Best_rime_rate
        Best_rime_rate=Rime_rates(1,i); Best_rime=Rimepop(i,:);
    end
end

T = 20;             % stagnation threshold
epsb = 1e-4;        % bucket for dedup

while FEs < MaxFEs
    RimeFactor = (rand-0.5)*2*cos((pi*FEs/(MaxFEs/10)))*(1-round(FEs*W/MaxFEs)/W);
    E = sqrt(FEs/MaxFEs);
    newRimepop = Rimepop;
    % normalized_rime_rates = normr(Rime_rates); % puncture disabled

    % Soft-rime only
    for i=1:size(Rimepop,1)
        for j=1:dim
            if rand()<E
                newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));
            end
        end
    end

    % Accept soft-rime moves greedily
    imp=false;
    for i=1:size(Rimepop,1)
        Flag4ub=newRimepop(i,:)>ub; Flag4lb=newRimepop(i,:)<lb;
        xi=(newRimepop(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        nf=fobj(xi); FEs=FEs+1;
        if nf<Rime_rates(1,i)
            Rime_rates(1,i)=nf; Rimepop(i,:)=xi; imp=true;
            if nf< Best_rime_rate, Best_rime_rate=nf; Best_rime=xi; end
        end
        if FEs>=MaxFEs, Convergence_curve(Time)=Best_rime_rate; return; end
    end

    % Update stagnation counter by global best
    if Best_rime_rate < lastBest
        noImp = 0; lastBest = Best_rime_rate;
    else
        noImp = noImp + 1;
    end

    % Occasional near-duplicate pruning (every few iters)
    if mod(Time, 10) == 0
        Q = floor(Rimepop/epsb);
        [~, first] = unique(Q, 'rows', 'stable');
        keep = false(size(Rimepop,1),1); keep(first) = true;
        Rimepop = Rimepop(keep,:); Rime_rates = Rime_rates(keep);
    end

    % GO miniblocks only if stagnated
    if noImp >= T && FEs < MaxFEs && size(Rimepop,1) >= 5
        [~, ind]=sort(Rime_rates,'ascend');
        % Adaptive elites and reflection prob with stagnation
        alpha = min(0.15, 0.10 + 0.005*(noImp - T));
        ke = max(2, round(alpha*size(Rimepop,1))); elite = ind(1:ke);
        topK = min(5, size(Rimepop,1));
        mut = 0.01 + (0.10-0.01)*(1 - FEs/MaxFEs);
        preflect = min(0.4, 0.2 + 0.01*(noImp - T));

        % GO learning on half elites
        L = elite(1:max(1,round(0.4*numel(elite))));
        Best = Rimepop(ind(1),:);
        for ii=1:numel(L)
            i=L(ii);
            Worst=Rimepop(ind(randi([max(1,size(Rimepop,1)-4), size(Rimepop,1)])),:);
            Better=Rimepop(ind(randi([2, min(5, size(Rimepop,1))])),:);
            v = [1:i-1, i+1:size(Rimepop,1)]; if numel(v) < 2, continue; end
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
                if rand<preflect
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

        % reset stagnation after a miniblock
        noImp = 0; lastBest = Best_rime_rate;
    end

    Convergence_curve(Time)=Best_rime_rate; Time=Time+1;
end
end

