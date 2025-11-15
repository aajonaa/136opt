function  Orderhao(filename)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    [~,~,rawdata] =  xlsread(filename,'overall');
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
    funcnum = [num(1) - 1]/(cur - 1);%%函数的数目
    %%函数的下标获取为 （k - 1）*（cur - 1） +2 —— （k - 1）*（cur - 1) + cur
    for i = 1:funcnum
        k = 1;
        for j = (i - 1) * (cur - 1) + 2 : (i - 1)*(cur - 1) + cur
            f(k) = cell2mat(rawdata(j,5));
            k = k + 1;
        end
        [~,index] = sort(f);
        level = 1;
        for j = 1:cur - 1
            fun(index(j)).func(i) = level ;
            level = level + 1;
            if j>1 && f(index(j)) == f(index(j-1))
                fun(index(j)).func(i) = fun(index(j-1)).func(i);
            end
        end
    end
    for i = 1 : cur - 1
        fun(i).mean = mean(fun(i).func);
        temp(i) = fun(i).mean;
    end
    [~,index] = sort(temp);
    result{1,1} = 'F';
    for i = 1:funcnum
        result{i+1,1} = cell2mat(rawdata((i - 1)*(cur - 1) +2,1));
    end
    result{funcnum+2, 1} = 'mean_level';
    result{funcnum+3, 1} = '平均结果';
    for i = 1: cur-1
        result{1,i+1} = fun(i).name;
        for j = 1:funcnum
            result{j+1,i+1} = fun(i).func(j);
        end
        result{funcnum+2,i+1} = (fun(i).mean);
        result{funcnum+3,index(i)+1} =i;
        if i>1 && temp(index(i)) == temp(index(i-1))
            result{funcnum+3,index(i)+1} = result{funcnum+3,index(i - 1)+1};
        end
    end
    xlswrite(filename,result,'result');
    
    % 合并pValue与result
    algorithm_num = cur - 1;
    xlswrite(filename,result(:,1),'result & pValue','A1');
    xlswrite(filename,result(:,2),'result & pValue','B1');
    for i = 1 : algorithm_num - 1
        
        number = 'B'+(i*3-2) - 'A';
        numberOfA = floor(number / 26);
        numberOfVar = mod(number, 26);

        frontName = [];
        if numberOfA ~= 0
            for j = 1 : numberOfA
                frontName = [frontName, 'A'];
            end
        end
        behindName = char('A' + numberOfVar);
        
        xlswrite(filename,result(:,i+2),'result & pValue',[[frontName, behindName],'1']);
    end
end

