function [Best_rime,Convergence_curve]=RIME_plusGOLearn(N,MaxFEs,lb,ub,dim,fobj)
% Alias naming to match RIME_plusGOReflect; wraps learning-only enhancement
[Best_rime,Convergence_curve]=RIME_plusLearn(N,MaxFEs,lb,ub,dim,fobj);
end

