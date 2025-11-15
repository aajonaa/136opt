function [best_pos, Convergence_curve]= GO_onlyReflect(N, MaxFEs, lb, ub, dim, fobj)
cfg = struct('do_learn',false,'do_reflect',true,'allow_random_accept',true,'allow_mutation',true);
[best_pos, Convergence_curve] = GO_core(N, MaxFEs, lb, ub, dim, fobj, cfg);
end

