function [best_pos,Convergence_curve]= GO_core(N, MaxFEs, lb, ub, dim, fobj, cfg)
% GO_core with ablation toggles
% cfg: do_learn(true), do_reflect(true), allow_random_accept(true), allow_mutation(true)
if nargin < 7 || isempty(cfg), cfg = struct(); end
cfg = set_default(cfg,'do_learn',true);
cfg = set_default(cfg,'do_reflect',true);
cfg = set_default(cfg,'allow_random_accept',true);
cfg = set_default(cfg,'allow_mutation',true);

FEs=0;
X=lb+(ub-lb)*unifrnd(0,1,N,dim);
Convergence_curve = [];
it = 1;
bestFitness=inf; AllFitness = zeros(1,N);
for i=1:N
    AllFitness(i)=fobj(X(i,:)); FEs=FEs+1;
    if bestFitness>AllFitness(i), bestFitness=AllFitness(i); best_pos=X(i,:); end
end
Convergence_curve(it)=bestFitness;
while FEs <= MaxFEs
    [~, ind]=sort(AllFitness);
    Best_X=X(ind(1),:);
    %% Learning phase
    if cfg.do_learn
        for i=1:N
            Worst_X = X(ind(randi([max(1,N-4),N])),:);
            Better_X=X(ind(randi([2,min(5,N)])),:);
            random=selectID(N,i,2);
            L1=random(1); L2=random(2);
            D1=(Best_X-Better_X); D2=(Best_X-Worst_X); D3=(Better_X-Worst_X); D4=(X(L1,:)-X(L2,:));
            s = norm(D1)+norm(D2)+norm(D3)+norm(D4);
            if s==0, s=eps; end
            LF1=norm(D1)/s; LF2=norm(D2)/s; LF3=norm(D3)/s; LF4=norm(D4)/s;
            SF=(AllFitness(i)/max(AllFitness));
            Gap=LF1*SF*D1 + LF2*SF*D2 + LF3*SF*D3 + LF4*SF*D4;
            newX=X(i,:)+Gap;
            newX=max(min(newX,ub),lb);
            newFit=fobj(newX); FEs=FEs+1;
            if AllFitness(i)>newFit || (cfg.allow_random_accept && rand<0.001 && ind(i)~=1)
                AllFitness(i)=newFit; X(i,:)=newX;
            end
            if bestFitness>AllFitness(i), bestFitness=AllFitness(i); best_pos=X(i,:); end
            if FEs>MaxFEs, break; end
        end
        if FEs>MaxFEs, break; end
    end

    %% Reflection phase
    if cfg.do_reflect
        for i=1:N
            newX=X(i,:);
            for j=1:dim
                if rand<0.3
                    R=X(ind(randi(min(5,N))),:);
                    newX(j)=X(i,j)+(R(:,j)-X(i,j))*unifrnd(0,1);
                    if cfg.allow_mutation && rand<(0.01+(0.1-0.01)*(1-FEs/MaxFEs))
                        newX(j)=lb+(ub-lb)*unifrnd(0,1);
                    end
                end
            end
            newX=max(min(newX,ub),lb);
            newFit=fobj(newX); FEs=FEs+1;
            if AllFitness(i)>newFit || (cfg.allow_random_accept && rand<0.001 && ind(i)~=1)
                AllFitness(i)=newFit; X(i,:)=newX;
            end
            if bestFitness>AllFitness(i), bestFitness=AllFitness(i); best_pos=X(i,:); end
            if FEs>MaxFEs, break; end
        end
        if FEs>MaxFEs, break; end
    end
    it = it + 1; Convergence_curve(it) = bestFitness;
end
end

function cfg = set_default(cfg, field, val)
if ~isfield(cfg,field) || isempty(cfg.(field)), cfg.(field)=val; end
end

%Subfunction:  Select two individuals different from Xi
function random = selectID(N, i, Num)
    if Num <= N
        vec=[1:i-1,i+1:N];
        r=zeros(1,2);
        for k =1:2
            n = N-k;
            t = randi(n,1,1);
            r(k) = vec(t);
            vec(t)=[];
        end
        random(1)=r(1);random(2)=r(2);
    else
        random = [max(1,mod(i,N)+1), max(1,mod(i+1,N)+1)];
    end
end

