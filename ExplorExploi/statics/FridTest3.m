function  FridTest3( xlsfilename ,fold2)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
readFilename=xlsfilename;     
writeFilename=xlsfilename;
[~,~,rawdata] = xlsread(xlsfilename,'overall');
num = size(rawdata);
flag = 0;
cur = 1;
for i = 2:num(1)
    if flag == 0
        fun(cur).name = cell2mat(rawdata(i,2));
        cur = cur + 1; 
    else
        break;
    end
    if i+1<=num(1) && ~strcmp(rawdata(i,1),rawdata(i+1,1))
        flag = 1;
    end
end

function_num = (num(1) - 1)/(cur - 1);         
algorithm_num = cur - 1;         
fold = fold2;                  

meanranks = [];
for i = 1: function_num
	dealData = xlsread(readFilename,['F',num2str(i)],'','basic');
	dealData = dealData(:,end); %%获取最后一列数据
    for j = 1:algorithm_num
	    data(:,j) = dealData((j-1)*fold+1:j*fold,:);
	end
	x = data;
	[p(i),table,stat(i)] = friedman(x);
	meanranks = [meanranks;stat(i).meanranks];
end


%%meanranks
meanranks = [meanranks;mean(meanranks)];%%排名
for i = 1:algorithm_num
	ans{1,i+1} = fun(i).name; 
    for j = 1:function_num+1
        ans{j+1,i+1} = meanranks(j,i);
    end
end
ans{1, algorithm_num + 2} = 'p';
ans{1,1} = 'F';

for i = 1:function_num
    ans{i+1,1} = cell2mat(rawdata((i - 1)*(cur - 1) +2,1));
    ans{i+1,algorithm_num + 2} = p(i);
end

ans{function_num+2, 1} = 'mean_level';
ans{function_num+3, 1} = '平均结果';
for i = 1:algorithm_num
    meanrank(i) = meanranks(function_num + 1,i);
end
[~,index] = sort(meanrank);
for i = 1:algorithm_num
   ans{function_num+3 ,index(i)+1} = i;
   if i>1 && meanrank(index(i)) == meanrank(index(i - 1))
       ans{function_num+3 ,index(i)+1} = ans{function_num+3 ,index(i-1)+1};
   end
end
xlswrite(xlsfilename,ans,'ans');

