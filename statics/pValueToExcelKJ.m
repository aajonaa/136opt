function  pValueToExcelKJ( xlsfilename ,fold2)
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
        cur = cur + 1; 
    else
        break;
    end
    if i+1<=num(1) && ~strcmp(rawdata(i,1),rawdata(i+1,1))
        flag = 1;
    end
end

functioon_num = (num(1) - 1)/(cur - 1);         
algorithm_num = cur - 1;         
fold = fold2;                  

sheet_name='pValue';
tag=['B3:B',num2str(1+algorithm_num)]; 
[~,label] = xlsread(readFilename,'overall',tag);
label=label';
xlswrite(writeFilename, label, sheet_name, 'B1')          %产生pValue表中的第一行

for i=1:functioon_num
    Function_name=['F',num2str(i)]; 
    Function_name={Function_name};
    xlswrite(writeFilename, Function_name, sheet_name, ['A',num2str(i+1)]);   %pValue表中第一列
end

for i=1:functioon_num 
    temp = xlsread(readFilename,['F',num2str(i)],'','basic');
    data(:,i) = temp(:,end);
end
pvalue=zeros(functioon_num,(algorithm_num-1));
for i=1:functioon_num
    for j=1:(algorithm_num-1)
        pvalue(i,j)=signrank(data(1:fold,i),data(j*fold+1:j*fold+fold,i));
    end
end
xlswrite(writeFilename, pvalue, sheet_name, 'B2')              %产生F1中的所有测试数据，从C1开始
end

