%% 协方差计算
function [ ZPosition ] = Covariance( RouBestPosition )
%COVARIANCE 此处显示有关此函数的摘要
%   此处显示详细说明
sigma=1.5;
W=zeros(size(RouBestPosition,1),size(RouBestPosition,2));
C=zeros(size(RouBestPosition,2),size(RouBestPosition,2));
M=mean(RouBestPosition);
for i=1:size(RouBestPosition,1)
    W(i,:)=RouBestPosition(i,:)-M;
end
for i=1:size(RouBestPosition,2)
    for j=1:size(RouBestPosition,2)
        C(i,j)=(sum(W(:,i).*W(:,j)))./(size(RouBestPosition,1)-1);
    end
end
[B,D]=eig(C);
r=rand(1,size(RouBestPosition,2));
ZPosition=M+(sigma*B*D*r')';
end

