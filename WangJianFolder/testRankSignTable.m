% xlsFileName = 'MSAO-0_26.xlsx';
% value = readcell(xlsFileName, 'sheet', 'pValue');
% algorithmName = {'MSAO', 'ALCPSO', 'BLPSO', 'CMAES', 'EOBLSSA', 'EPSDE',...
%     'GOTLBO', 'CBA', 'CLPSO', 'CAGWO', 'GWOCMA', 'RUN', 'CCMWOA', 'RDWOA', ...
%     'SMA', 'RIME', 'HGS', 'INFO'} %#ok<NOPTS>
% funcNum = 12;
% algoNum = size(algorithmName, 2);
%% Write the pValue & sign table
pValueSign{1, 1} = 'Func';
for func = 1:funcNum % The first column
    pValueSign{func+1, 1} = ['F', num2str(func)];
end
for algo = 1:algoNum-1 % The first row
    pValueSign{1, algo+1} = algorithmName{algo+1};
end
% RankPValue = xlsread(xlsFileName, 'rank & pValue');
RankPValue = readmatrix(xlsFileName, 'sheet', 'rank & pValue');
RankPValue = RankPValue(:, 2:end);
pValueSign = cell(funcNum+2, algoNum);
pValueSign{1, 1} = 'Func';
pValueSign{funcNum+2, 1} = '+/-/=';
pValueSign{funcNum+3, 1} = '+/-/=';
for func = 1:funcNum % The first column
    pValueSign{func+1, 1} = num2str(func);
end
for algo = 1:algoNum-1 % The first row
    pValueSign{1, 2*algo} = algorithmName{algo+1};
end
for algo = 1:algoNum-1 % The pValue and sign table
    cntBetter = 0;
    cntWorse = 0;
    for func = 1:funcNum
        if RankPValue(func, 2*algo+1) < 0.05
            if RankPValue(func, 1) < RankPValue(func, 2*algo)
                pValueSign{func+1, 2*algo+1} = '+';
                cntBetter = cntBetter + 1;
            else
                pValueSign{func+1, 2*algo+1} = '-';
                cntWorse = cntWorse + 1;
            end
        else
            pValueSign{func+1, 2*algo+1} = '=';
        end
        pValueSign{func+1, 2*algo} = num2str(value{func+1, algo+1});
    end
    cntEqual = funcNum - cntBetter - cntWorse;
    pValueSign{func+2, 2*algo} = [num2str(cntBetter), '/', num2str(cntWorse), ...
        '/', num2str(cntEqual)];
    pValueSign{func+3, algo+1} = [num2str(cntBetter), '/', num2str(cntWorse), ...
        '/', num2str(cntEqual)];
end
writecell(pValueSign, xlsFileName, 'sheet', 'pValue & sign');