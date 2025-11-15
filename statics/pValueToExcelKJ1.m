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
pvalue=zeros(algorithm_num-1, algorithm_num, functioon_num);
for i = 1 : algorithm_num-1
    for j = i+1 : algorithm_num
        for k = 1 : functioon_num
            pvalue(i,j,k)=signrank(data((i-1)*fold+1 : (i-1)*fold+fold, k), data((j-1)*fold+1 : (j-1)*fold+fold, k));
        end
    end
end
for i = 1 : algorithm_num-1
    pvalueW = zeros(functioon_num, algorithm_num-i);
    for j = i+1 : algorithm_num
        pvalueW(:, j-i) = pvalue(i,j,:);
    end
    xlswrite(writeFilename, pvalueW, sheet_name, [char('B' + (i-1)), num2str((i-1)*functioon_num+2)])              %产生F1中的所有测试数据，从C1开始
end

% 合并pValue与result
tmpdata = xlsread(writeFilename, 'pValue');
for i = 1 : algorithm_num - 1
    for j = i+1 : algorithm_num
    
        number = 'C'+((j-1)*3-2) - 'A';
        numberOfA = floor(number / 26);
        numberOfVar = mod(number, 26);

        frontName = [];
        if numberOfA ~= 0
            for k = 1 : numberOfA
                frontName = [frontName, 'A'];
            end
        end
        behindName = char('A' + numberOfVar);

        xlswrite(writeFilename, tmpdata((i-1)*functioon_num+1:(i-1)*functioon_num+functioon_num, j-1), 'result & pValue',[[frontName, behindName], num2str((i-1)*(functioon_num+5) + 2)]);
        xlswrite(writeFilename, {'pvalue'}, 'result & pValue',[[frontName, behindName],'1']);
    
    end
end
% 统计强弱
xlswrite(writeFilename, {'CountOfBetter'}, 'result & pValue',['A',num2str(functioon_num + 4)]);
xlswrite(writeFilename, {'CountOfWorse'}, 'result & pValue',['A',num2str(functioon_num + 5)]);
xlswrite(writeFilename, {'CountOfEqual'}, 'result & pValue',['A',num2str(functioon_num + 6)]);
tmpAlpha = ' ';
tmpResult = repmat(tmpAlpha, (functioon_num+5)*(algorithm_num - 1), algorithm_num - 1);
tmpdata = xlsread(writeFilename, 'result & pValue');
for i = 1 : algorithm_num - 1
    for k = i+1 : algorithm_num
        
        cntOfBetter = 0;
        cntOfWorse = 0;
        for j = 1 : functioon_num
            if tmpdata((i-1)*functioon_num+j, (k-1)*3) < 0.05
                if (i==1)
                    tmpIndex = 1;
                else
                    tmpIndex = i*3-4;
                end
                if tmpdata(j, tmpIndex) < tmpdata(j, (k-1)*3-1)
                    tmpResult((i-1)*(functioon_num+5)+j, k-1) = '+';
                    cntOfBetter = cntOfBetter + 1;
                else
                    tmpResult((i-1)*(functioon_num+5)+j, k-1) = '-';
                    cntOfWorse = cntOfWorse + 1;
                end
            end
        end
        cntOfEqual = functioon_num - cntOfBetter - cntOfWorse;

        number = 'C'+((k-1)*3-1) - 'A';
        numberOfA = floor(number / 26);
        numberOfVar = mod(number, 26);

        frontName = [];
        if numberOfA ~= 0
            for j = 1 : numberOfA
                frontName = [frontName, 'A'];
            end
        end
        behindName = char('A' + numberOfVar);


        xlswrite(writeFilename, tmpResult((i-1)*(functioon_num+5)+1:(i-1)*(functioon_num+5)+functioon_num, (k-1)), 'result & pValue',[[frontName, behindName], num2str((i-1)*(functioon_num+5) + 2)]);
        xlswrite(writeFilename, cntOfBetter, 'result & pValue',[[frontName, behindName],num2str((i-1)*(functioon_num+5) + functioon_num + 4)]);
        xlswrite(writeFilename, cntOfWorse, 'result & pValue',[[frontName, behindName],num2str((i-1)*(functioon_num+5) + functioon_num + 5)]);
        xlswrite(writeFilename, cntOfEqual, 'result & pValue',[[frontName, behindName],num2str((i-1)*(functioon_num+5) + functioon_num + 6)]);
        
    end
end

