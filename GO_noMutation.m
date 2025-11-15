function [best_pos, Convergence_curve]= GO_noMutation(N, MaxFEs, lb, ub, dim, fobj)
cfg = struct('do_learn',true,'do_reflect',true,'allow_random_accept',true,'allow_mutation',false);
[best_pos, Convergence_curve] = GO_core(N, MaxFEs, lb, ub, dim, fobj, cfg);
end

