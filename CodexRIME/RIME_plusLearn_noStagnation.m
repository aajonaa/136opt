function [Best_rime,Convergence_curve]=RIME_plusLearn_noStagnation(N,MaxFEs,lb,ub,dim,fobj)
cfg = struct('do_soft',true,'do_punct',false,'stagnation_gate',false, ...
             'T',20,'alpha',0.10,'learnFrac',0.40,'topK',5, ...
             'includeD4',true,'useSF',true);
[Best_rime,Convergence_curve]=RIME_plusLearn_core(N,MaxFEs,lb,ub,dim,fobj,cfg);
end

