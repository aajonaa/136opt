function [best_pos, Convergence_curve]= GO_noRandomAccept(N, MaxFEs, lb, ub, dim, fobj)
cfg = struct('do_learn',true,'do_reflect',true,'allow_random_accept',false,'allow_mutation',true);
[best_pos, Convergence_curve] = GO_core(N, MaxFEs, lb, ub, dim, fobj, cfg);
end

