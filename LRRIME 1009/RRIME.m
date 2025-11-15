function [Best_rime,Convergence_curve]=RIME_plusGOReflect(N,MaxFEs,lb,ub,dim,fobj)
% RIME with GO-style reflection miniblock, guided by ablation:
% - Keep RIME soft-rime (dominant)
% - Disable hard puncture (less helpful)
% - Add GO reflection on elites
% - No random acceptance; keep mutation (decaying)

Best_rime=zeros(1,dim);
Best_rime_rate=inf;
Rimepop=initialization(N,dim,ub,lb);
Lb=lb.*ones(1,dim); Ub=ub.*ones(1,dim);
FEs=0; Time = 1; Convergence_curve=[];
Rime_rates=zeros(1,N); newRime_rates=zeros(1,N);
W = 5; noImp = 0;

% init eval
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

    % Soft-rime only (puncture disabled per ablation)
    for i=1:N
        for j=1:dim
            if rand()<E
                newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));
            end
            % no hard puncture
            % if rand()<normalized_rime_rates(i), newRimepop(i,j)=Best_rime(1,j); end
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

    % GO reflection miniblock on elites (no random acceptance)
    [~, ind]=sort(Rime_rates,'ascend');
    ke = max(2, round(0.10*N)); elite = ind(1:ke);
    topK = min(5,N);
    mut = 0.01 + (0.10-0.01)*(1 - FEs/MaxFEs);

    % Trigger each iter or when stagnated; cheap and bounded
    if (true || noImp>=20) && FEs<MaxFEs
        for ii=1:numel(elite)
            i = elite(ii);
            cand = Rimepop(i,:);
            for j=1:dim
                if rand<0.30
                    R = Rimepop(ind(randi(topK)),:);
                    cand(j)=Rimepop(i,j) + (R(j)-Rimepop(i,j))*rand;
                    if rand<mut, cand(j)=lb + (ub-lb).*rand(1,1); end
                end
            end
            % clip, eval, greedy accept only
            cand = max(min(cand,ub),lb);
            if FEs>=MaxFEs, break; end
            nf = fobj(cand); FEs=FEs+1;
            if nf < Rime_rates(1,i)
                Rime_rates(1,i) = nf; Rimepop(i,:) = cand; imp=true;
                if nf < Best_rime_rate, Best_rime_rate = nf; Best_rime = cand; end
            end
            if FEs>=MaxFEs, break; end
        end
    end

    Convergence_curve(Time)=Best_rime_rate; Time=Time+1;
end
end

