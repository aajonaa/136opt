function [best_pos, curve] = RIME_GO_core(N, MaxFEs, lb, ub, dim, fobj, cfg)
% RIME_GO_core - Hybrid RIME + GO with configurable ablations via cfg
% cfg fields (all logical unless noted):
%   do_rime_step (true)
%   do_punct (true)
%   do_coeff_move (true)
%   do_prune_dup (true)
%   do_prune_agefit (true)
%   do_adaptive_N (true)
%   do_GO (true)
%   do_GO_learn (true)
%   do_GO_reflect (true)

if nargin < 7 || isempty(cfg)
    cfg = struct();
end
cfg = set_default(cfg, 'do_rime_step', true);
cfg = set_default(cfg, 'do_punct', true);
cfg = set_default(cfg, 'do_coeff_move', true);
cfg = set_default(cfg, 'do_prune_dup', true);
cfg = set_default(cfg, 'do_prune_agefit', true);
cfg = set_default(cfg, 'do_adaptive_N', true);
cfg = set_default(cfg, 'do_GO', true);
cfg = set_default(cfg, 'do_GO_learn', true);
cfg = set_default(cfg, 'do_GO_reflect', true);

if isscalar(lb), lb = lb.*ones(1,dim); end
if isscalar(ub), ub = ub.*ones(1,dim); end

% init
X = repmat(lb,N,1) + rand(N,dim).*repmat((ub-lb),N,1);
F = inf(1,N); ages = zeros(1,N); W = 5; it=1; curve=[]; FEs=0;
bestFitness = inf; best_pos = zeros(1,dim); N0=N;

% eval init
for i=1:N
    F(i)=fobj(X(i,:)); FEs=FEs+1;
    if F(i)<bestFitness, bestFitness=F(i); best_pos=X(i,:); end
    if FEs>=MaxFEs, curve(it)=bestFitness; return; end
end

noImp=0;
while FEs<MaxFEs
    % --- RIME step ---
    coeff = (rand-0.5)*2*cos((pi*FEs/(MaxFEs/10)))*(1-round(FEs*W/MaxFEs)/W);
    p = sqrt(FEs/MaxFEs);
    fmin=min(F); fmax=max(F); den=(fmax-fmin)+eps;
    punct = (F - fmin)/den;
    newX = X;
    if cfg.do_rime_step
        for i=1:size(X,1)
            for j=1:dim
                if cfg.do_coeff_move && rand<p
                    newX(i,j)=best_pos(j)+coeff*((ub(j)-lb(j))*rand+lb(j));
                end
                if cfg.do_punct && rand<punct(i)
                    newX(i,j)=best_pos(j);
                end
            end
        end
    end

    % accept
    imp=false;
    for i=1:size(X,1)
        xi = min(max(newX(i,:),lb),ub);
        nf = fobj(xi); FEs=FEs+1;
        if nf<F(i)
            X(i,:)=xi; F(i)=nf; ages(i)=0; imp=true;
            if nf<bestFitness, bestFitness=nf; best_pos=xi; end
        else
            ages(i)=ages(i)+1;
        end
        if FEs>=MaxFEs, curve(it)=bestFitness; return; end
    end
    noImp = imp*0 + (~imp)*(noImp+1);

    % --- PRUNE (near-dup + age×fitness) ---
    if cfg.do_prune_dup
        epsb = 1e-4; Q = floor(X/epsb);
        [~, first] = unique(Q, 'rows', 'stable');
        keep = false(size(X,1), 1);
        keep(first) = true;
        X = X(keep, :); F = F(keep); ages = ages(keep);
    end

    if cfg.do_prune_agefit
        [~,rf]=sort(F,'ascend'); rF=zeros(size(rf)); rF(rf)=1:numel(rf);
        [~,ra]=sort(ages,'ascend'); rA=zeros(size(ra)); rA(ra)=1:numel(ra);
        score = rF + 0.6*rA;
        frac = 0.10 + 0.20*(FEs/MaxFEs);
        K = max(0, round(frac*size(X,1)));
        if K>0
            [~,ord]=sort(score,'ascend');
            idx=ord(1:end-K);
            X=X(idx,:); F=F(idx); ages=ages(idx);
        end
        X = min(max(X,lb),ub);
    end

    % --- Adaptive N + refill ---
    if cfg.do_adaptive_N
        Nend = max(5, round(0.60*N0));
        Ntgt = max(5, round(N0 - (N0 - Nend)*(FEs/MaxFEs)));
        need = Ntgt - size(X,1);
        if need>0 && FEs<MaxFEs
            Xr = refill(need,dim,lb,ub);
            Fr = zeros(1,size(Xr,1));
            kkeep = false(1,size(Xr,1));
            for r=1:size(Xr,1)
                if FEs>=MaxFEs, break; end
                Fr(r)=fobj(Xr(r,:)); FEs=FEs+1; kkeep(r)=true;
                if Fr(r)<bestFitness, bestFitness=Fr(r); best_pos=Xr(r,:); end
            end
            if any(kkeep)
                X=[X; Xr(kkeep,:)]; F=[F, Fr(kkeep)]; ages=[ages, zeros(1,sum(kkeep))];
            end
        end
    end

    % --- Tiny GO on elites ---
    if cfg.do_GO && (mod(it,1)==0 || noImp>=30) && size(X,1)>=5 && FEs<MaxFEs
        [F,X,ages,bestFitness,best_pos,FEs] = go_miniblock(F,X,ages,lb,ub,bestFitness,best_pos,fobj,FEs,MaxFEs,cfg);
        noImp=0;
    end

    curve(it)=bestFitness; it=it+1;
end

end

% ======== helpers ========
function cfg = set_default(cfg, field, val)
if ~isfield(cfg, field) || isempty(cfg.(field))
    cfg.(field) = val;
end
end

function Xr = refill(m,dim,lb,ub)
% Memory-safe space-filling refill using Latin Hypercube Sampling (LHS)
if m<=0, Xr=zeros(0,dim); return; end
U = zeros(m, dim);
for j = 1:dim
    perm = randperm(m).';
    U(:, j) = (perm - rand(m,1)) / m;
end
U = min(max(U, 0), 1);
Xr = repmat(lb, m, 1) + U .* repmat((ub - lb), m, 1);
end

function [F,X,ages,bfit,bpos,FEs] = go_miniblock(F,X,ages,lb,ub,bfit,bpos,fobj,FEs,MaxFEs,cfg)
N=size(X,1); dim=size(X,2); [~,ind]=sort(F,'ascend');
ke = max(2, round(0.10*N)); elite = ind(1:ke);
topK = min(5,N);
mut = 0.01 + (0.10-0.01)*(1 - FEs/MaxFEs);

% learning (half elites)
if cfg.do_GO_learn
    L = elite(1:max(1,round(0.4*numel(elite))));
    Best = X(ind(1),:);
    for ii=1:numel(L)
        i=L(ii);
        Worst=X(ind(randi([max(1,N-4),N])),:);
        Better=X(ind(randi([2,min(5,N)])),:);
        [l1,l2]=pick2(N,i);
        D1=Best-Better; D2=Best-Worst; D3=Better-Worst; D4=X(l1,:)-X(l2,:);
        s = norm(D1)+norm(D2)+norm(D3)+norm(D4)+eps;
        Gap = (norm(D1)/s)*D1 + (norm(D2)/s)*D2 + (norm(D3)/s)*D3 + (norm(D4)/s)*D4;
        cand = min(max(X(i,:)+ (F(i)/(max(F)+eps))*Gap, lb), ub);
        if FEs>=MaxFEs, return; end
        nf=fobj(cand); FEs=FEs+1;
        if nf<F(i), X(i,:)=cand; F(i)=nf; ages(i)=0; if nf<bfit, bfit=nf; bpos=cand; end
        else, if rand<0.001 && i~=ind(1), X(i,:)=cand; F(i)=nf; ages(i)=ages(i)+1; if nf<bfit, bfit=nf; bpos=cand; end
            else, ages(i)=ages(i)+1; end
        end
    end
end

% reflection (other half)
if cfg.do_GO_reflect
    R = elite(max(1, numel(elite)-max(1,round(0.4*numel(elite)))+1):end);
    for ii=1:numel(R)
        i=R(ii); cand=X(i,:);
        for j=1:dim
            if rand<0.30
                Rv=X(ind(randi(topK)),:);
                cand(j)=X(i,j)+(Rv(j)-X(i,j))*rand;
                if rand<mut, cand(j)=lb(j)+(ub(j)-lb(j))*rand; end
            end
        end
        cand=min(max(cand,lb),ub);
        if FEs>=MaxFEs, return; end
        nf=fobj(cand); FEs=FEs+1;
        if nf<F(i), X(i,:)=cand; F(i)=nf; ages(i)=0; if nf<bfit, bfit=nf; bpos=cand; end
        else, if rand<0.001 && i~=ind(1), X(i,:)=cand; F(i)=nf; ages(i)=ages(i)+1; if nf<bfit, bfit=nf; bpos=cand; end
            else, ages(i)=ages(i)+1; end
        end
    end
end
end

function [a,b]=pick2(N,excl)
v=[1:excl-1, excl+1:N]; r1=v(randi(numel(v))); v(v==r1)=[];
r2=v(randi(numel(v))); a=r1; b=r2;
end

