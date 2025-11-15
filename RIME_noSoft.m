function [best_pos, curve] = RIME_noSoft(N,MaxFEs,lb,ub,dim,fobj)
cfg = struct('do_soft', false, 'do_punct', true);
[best_pos, curve] = RIME_core(N,MaxFEs,lb,ub,dim,fobj,cfg);
end

