function [ Total_divergence ] = Divergence_calculation( X,dim,N )
Divergence=zeros(dim,1);
Dis_sum=0;
for j=1:dim
    X(:,j)=sort(X(:,j));
    median=X(round(N/2),j);
    for i=1:N
        Dis_sum=Dis_sum+abs(median-X(i,j));
    end
    Divergence(j)=Dis_sum/N;
end
Total_divergence=sum(Divergence)/dim;
end

