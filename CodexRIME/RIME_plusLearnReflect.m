function [Best_rime,Convergence_curve]=RIME_plusLearnReflect(N,MaxFEs,lb,ub,dim,fobj)
% Alias wrapper for combined learning+reflection variant
[Best_rime,Convergence_curve]=RIME_plusGOLearnReflect(N,MaxFEs,lb,ub,dim,fobj);
end

