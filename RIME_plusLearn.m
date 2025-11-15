function [Best_rime,Convergence_curve]=RIME_plusLearn(N,MaxFEs,lb,ub,dim,fobj)
% Default: soft-rime on; puncture off; stagnation gated; elite 10%, learn 40%, topK=5
cfg = struct('do_soft',true,'do_punct',false,'stagnation_gate',true, ...
             'T',20,'alpha',0.10,'learnFrac',0.40,'topK',5, ...
             'includeD4',true,'useSF',true);
[Best_rime,Convergence_curve]=RIME_plusLearn_core(N,MaxFEs,lb,ub,dim,fobj,cfg);
end

