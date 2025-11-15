function [ L,M ] = Get_QDOA( Q,D )
%输入：水平数Q，维数D。
%输出：正交表L，正交表行数M。
J=floor(log(D*(Q-1))/log(Q)+1);
N=(Q.^J-1)/(Q-1);
M=Q.^J;
L=zeros(M,N);
for k=1:J
    j=(Q.^(k-1)-1)/(Q-1)+1;
    for i=1:M
        L(i,j)=mod(floor((i-1)/(Q.^(J-k))),Q);
    end
end
for k=2:J
    j=(Q.^(k-1)-1)/(Q-1)+1;
    for s=1:(j-1)
        for t=1:(Q-1)
            for i=1:M
                L(i,j+t+(s-1)*(Q-1))=mod(L(i,s)*t+L(i,j),Q);
            end
        end
    end
end
L=L(:,1:D);
L=L+1;
end

