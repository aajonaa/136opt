function [best_pos, curve] = RIME_GO_onlyReflect(N, MaxFEs, lb, ub, dim, fobj)
cfg = struct('do_rime_step',true,'do_punct',true,'do_coeff_move',true, ...
             'do_prune_dup',true,'do_prune_agefit',true, ...
             'do_adaptive_N',true,'do_GO',true,'do_GO_learn',false,'do_GO_reflect',true);
[best_pos, curve] = RIME_GO_core(N, MaxFEs, lb, ub, dim, fobj, cfg);
end

