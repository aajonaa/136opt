function [ k1,k2,k3 ] = GetRan3( num,N )
%从1到N中随机产生2个不等于num的整数,赋值给k1和k2
%   此处显示详细说明
k1=randperm(N,1);
while(k1==num)
  k1=randperm(N,1);
end
k2=randperm(N,1);
while(k2==num || k2==k1)
  k2=randperm(N,1);
end
k3=randperm(N,1);
while(k3==num || k2==k1 || k3==k2)
  k3=randperm(N,1);
end
end

