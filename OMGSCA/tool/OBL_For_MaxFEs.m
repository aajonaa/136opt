%% 反向学习
function [ P , fes ] = OBL_For_MaxFEs( X,N,dim,ub,lb,fobj,type,fes)
%    X种群，N种群大小，dim问题维度，lb下界，ub上界，fobj适应度函数，type反向学习模式选择参数，0为固定边界，1位动态边界。
%   默认上下边界向量是行方向，如果实际情况是列方向，请将上下边界向量转置后使用
ub=ones(1,dim).*ub;%转换上界向量为1行dim列向量
lb=ones(1,dim).*lb;%转换下界向量为1行dim列向量
S=zeros(N,dim);
if type==0 %%一般随机反向学习
    for i = 1:N
         k=rand();
         for j = 1:dim
            S(i,j) = k * (ub(j)+lb(j)) - X(i,j);
            if S(i,j) < lb(j)||S(i,j) > ub(j)
                S(i,j) = lb(j) + rand() * (ub(j) - lb(j));  
            end
         end
    end
	Positions1 = S;
	Positions2 = [X; Positions1];
	for i = 1:2*N
		fit_value(i) = fobj(Positions2(i, :));
        fes = fes + 1;
	end
	[fit_value,index] = sort(fit_value);
	P = Positions2(index, :); 
	P = P(1:N, :); 
elseif type==1 %%动态边界随机反向学习
    dub=max(X);
    dlb=min(X);
     for i = 1:N
         k=rand();
         for j = 1:dim
            S(i,j) = k * (dub(j)+dlb(j)) - X(i,j);
            if S(i,j) < lb(j)||S(i,j) > ub(j)
                S(i,j) = lb(j) + rand() * (ub(j) - lb(j));  
            end
         end 
     end 
    Positions1 = S;
    Positions2 = [X; Positions1];
    for i = 1:2*N
        fit_value(i) = fobj(Positions2(i, :));
        fes = fes + 1;
    end
    [fit_value,index] = sort(fit_value);
    P = Positions2(index, :); 
    P = P(1:N, :); 
end

end

