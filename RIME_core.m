function [Best_rime,Convergence_curve]=RIME_core(N,MaxFEs,lb,ub,dim,fobj,cfg)
% RIME_core with ablation toggles in cfg
% cfg fields: do_soft (true), do_punct (true)
if nargin < 7 || isempty(cfg), cfg = struct(); end
if ~isfield(cfg,'do_soft') || isempty(cfg.do_soft), cfg.do_soft = true; end
if ~isfield(cfg,'do_punct') || isempty(cfg.do_punct), cfg.do_punct = true; end

Best_rime=zeros(1,dim);
Best_rime_rate=inf;
Rimepop=initialization(N,dim,ub,lb);
Lb=lb.*ones(1,dim);
Ub=ub.*ones(1,dim);
FEs=0; Time = 1; Convergence_curve=[];
Rime_rates=zeros(1,N); newRime_rates=zeros(1,N);
W = 5;

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
    normalized_rime_rates = normr(Rime_rates); % requires NN toolbox; present in repo

    for i=1:N
        for j=1:dim
            if cfg.do_soft && rand()<E
                newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));
            end
            if cfg.do_punct && rand()<normalized_rime_rates(i)
                newRimepop(i,j)=Best_rime(1,j);
            end
        end
    end

    for i=1:N
        Flag4ub=newRimepop(i,:)>ub; Flag4lb=newRimepop(i,:)<lb;
        newRimepop(i,:)=(newRimepop(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        newRime_rates(1,i)=fobj(newRimepop(i,:)); FEs=FEs+1;
        if newRime_rates(1,i)<Rime_rates(1,i)
            Rime_rates(1,i) = newRime_rates(1,i);
            Rimepop(i,:) = newRimepop(i,:);
            if newRime_rates(1,i)< Best_rime_rate
               Best_rime_rate=Rime_rates(1,i); Best_rime=Rimepop(i,:);
            end
        end
    end
    Convergence_curve(Time)=Best_rime_rate; Time=Time+1;
end
end

