function [ out_ver ] = MySampling(in_ver,num)
in_num=size(in_ver,2);
out_ver=zeros(1,num);
step = round(in_num/num);
if step>0  % 输入数组比输出数组长
    out_ver=in_ver(round(linspace(1,in_num,num)));
else      % 输出数组比输入数组长，插入数据
    index=1:1:num;
    newLocalIndex=round(linspace(1,num,in_num));
    restIndex=setdiff(index, newLocalIndex);
    out_ver(newLocalIndex)=in_ver;
    for i=1:size(restIndex,2)
       if restIndex(i) == 1
           restIndex(i) = 2;
       end
       out_ver(restIndex(i))=out_ver(restIndex(i)-1);
    end
%     out_ver(restIndex)=out_ver(restIndex-1);

end
end

