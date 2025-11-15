function [best_pos, curve] = RIME_noPunct(N,MaxFEs,lb,ub,dim,fobj)
cfg = struct('do_soft', true, 'do_punct', false);
[best_pos, curve] = RIME_core(N,MaxFEs,lb,ub,dim,fobj,cfg);
end

