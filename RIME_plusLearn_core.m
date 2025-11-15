function [Best_rime,Convergence_curve]=RIME_plusLearn_core(N,MaxFEs,lb,ub,dim,fobj,cfg)
% RIME with GO-style learning miniblock; configurable for ablation
% cfg fields (defaults in parentheses):
%   do_soft (true)        - RIME soft-rime step
%   do_punct (false)      - RIME puncture (recommended false per ablation)
%   stagnation_gate (true)- gate learning on stagnation
%   T (20)                - stagnation threshold (iters without global improv.)
%   alpha (0.10)          - elite fraction for learning
%   learnFrac (0.40)      - fraction of elites to apply learning
%   topK (5)              - top-K bound for "Better"/"Worst" picks
%   includeD4 (true)      - include random pair difference term D4
%   useSF (true)          - use SF fitness scaling on the gap

if nargin < 7 || isempty(cfg), cfg = struct(); end
cfg = setdef(cfg,'do_soft',true);
cfg = setdef(cfg,'do_punct',false);
cfg = setdef(cfg,'stagnation_gate',true);
cfg = setdef(cfg,'T',20);
cfg = setdef(cfg,'alpha',0.10);
cfg = setdef(cfg,'learnFrac',0.40);
cfg = setdef(cfg,'topK',5);
cfg = setdef(cfg,'includeD4',true);
cfg = setdef(cfg,'useSF',true);

Best_rime=zeros(1,dim); Best_rime_rate=inf;
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

while FEs < MaxFEs
    RimeFactor = (rand-0.5)*2*cos((pi*FEs/(MaxFEs/10)))*(1-round(FEs*W/MaxFEs)/W);
    E = sqrt(FEs/MaxFEs);
    newRimepop = Rimepop;
    if cfg.do_punct
        normalized_rime_rates = normr(Rime_rates); %#ok<UNRCH>
    end

    % RIME soft/punct steps
    for i=1:size(Rimepop,1)
        for j=1:dim
            if cfg.do_soft && rand()<E
                newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));
            end
            if cfg.do_punct && rand()<normalized_rime_rates(i)
                newRimepop(i,j)=Best_rime(1,j);
            end
        end
    end

    % Accept RIME moves greedily
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

    % Update stagnation
    if Best_rime_rate < lastBest
        noImp = 0; lastBest = Best_rime_rate;
    else
        noImp = noImp + 1;
    end

    % Learning miniblock: gated or always
    if (~cfg.stagnation_gate || noImp >= cfg.T) && FEs < MaxFEs && size(Rimepop,1) >= 5
        [~, ind]=sort(Rime_rates,'ascend');
        ke = max(2, round(cfg.alpha*size(Rimepop,1))); elite = ind(1:ke);
        L = elite(1:max(1,round(cfg.learnFrac*numel(elite))));
        Best = Rimepop(ind(1),:);
        topK = min(cfg.topK, size(Rimepop,1));

        for ii=1:numel(L)
            i=L(ii);
            Worst=Rimepop(ind(randi([max(1,size(Rimepop,1)-4), size(Rimepop,1)])),:);
            Better=Rimepop(ind(randi([2, topK])),:);
            % random distinct pair for D4 if enabled
            if cfg.includeD4
                v = [1:i-1, i+1:size(Rimepop,1)]; if numel(v) < 2, continue; end
                r1=v(randi(numel(v))); v(v==r1)=[]; r2=v(randi(numel(v)));
                D4=Rimepop(r1,:)-Rimepop(r2,:);
            else
                D4=zeros(1,dim);
            end
            D1=Best-Better; D2=Best-Worst; D3=Better-Worst;
            s = norm(D1)+norm(D2)+norm(D3)+norm(D4); if s==0, s=eps; end
            % weights proportional to contribution norms
            Gap = (norm(D1)/s)*D1 + (norm(D2)/s)*D2 + (norm(D3)/s)*D3 + (norm(D4)/s)*D4;
            if cfg.useSF
                SF = (Rime_rates(i)/max(Rime_rates));
                cand = Rimepop(i,:) + SF * Gap;
            else
                cand = Rimepop(i,:) + Gap;
            end
            cand = max(min(cand,ub),lb);
            if FEs>=MaxFEs, break; end
            nf=fobj(cand); FEs=FEs+1;
            if nf<Rime_rates(1,i)
                Rime_rates(1,i)=nf; Rimepop(i,:)=cand; imp=true;
                if nf< Best_rime_rate, Best_rime_rate=nf; Best_rime=cand; end
            end
            if FEs>=MaxFEs, break; end
        end

        % reset stagnation if gated
        if cfg.stagnation_gate, noImp = 0; lastBest = Best_rime_rate; end
    end

    Convergence_curve(Time)=Best_rime_rate; Time=Time+1;
end
end

function s = setdef(s, k, v)
if ~isfield(s,k) || isempty(s.(k)), s.(k) = v; end
end

