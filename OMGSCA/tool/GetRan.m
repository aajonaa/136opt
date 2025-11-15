function [ k ] = GetRan( num,N )
%随机产生一个不等于num的整数
%   此处显示详细说明
k=randperm(N,1);
while(k==num)
  k=randperm(N,1);
end
end

