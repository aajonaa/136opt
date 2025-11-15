function [best_pos, curve] = RIME_GO(N, MaxFEs, lb, ub, dim, fobj)
% Thin wrapper calling the configurable core with default settings
cfg = [];
[best_pos, curve] = RIME_GO_core(N, MaxFEs, lb, ub, dim, fobj, cfg);
end
