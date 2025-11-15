 % for i = 1:20
%% V2 version add the box plot
%% 2025-2-11 代码存在某些问题：AvgFitnessRank:
% SBO : 1
% AvgFitnessRankSum:
% SBO : 1
% AvgRank:
% SBO : 0.034483
% Rank:
% SBO : 1
%% Initial setting
clear; close all;  clc; tic;
%% 1.Legend off
%% 2.Change FunctionSeq name

%% 
% 将仓库根目录与脚本目录加入路径，确保能找到 Get_Functions 等方法
% 注意：脚本中 mfilename('fullpath') 可能返回空，需要用 dbstack 作为兜底
stackInfo = dbstack('-completenames');
if ~isempty(stackInfo)
    scriptDir = fileparts(stackInfo(1).file);
else
    % 作为保底，使用 which；若仍为空，则使用当前工作目录
    scriptPathGuess = which('ImproveTest_v2');
    if ~isempty(scriptPathGuess)
        scriptDir = fileparts(scriptPathGuess);
    else
        scriptDir = pwd;
    end
end
projectRoot = fileparts(scriptDir);
addpath(genpath(projectRoot));
addpath(genpath(scriptDir));

%% Test functions choose
% 'CEC2017' %#ok<NOPTS>
% FunctionSeq is the FunctionTittle of the figure
% FunctionSeq       = {'F1',  'F3',  'F4',  'F5',  'F6',  'F7',  'F8',  'F9',  'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18', 'F19', 'F20', 'F21', 'F22', 'F23', 'F24', 'F25', 'F26', 'F27', 'F28', 'F29', 'F30'};
FunctionSeq       = {'F1',  'F2',  'F3',  'F4',  'F5',  'F6',  'F7',  'F8',  'F9',  'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18', 'F19', 'F20', 'F21', 'F22', 'F23', 'F24', 'F25', 'F26', 'F27', 'F28', 'F29'};
Function_name_all = {'F107','F109','F110','F111','F112','F113','F114','F115','F116','F117','F118','F119','F120','F121','F122','F123','F124','F125','F126','F127','F128','F129','F130','F131','F132','F133','F134','F135','F136'};
% FunctionSeq       = {'F26'};
% Function_name_all = {'F133'};
% 'CEC2019' %#ok<NOPTS>
% FunctionSeq       = {'F1',  'F2',  'F3',  'F4',  'F5',  'F6',  'F7',  'F8',  'F9',  'F10'};
% Function_name_all = {'F137','F138','F139','F140','F141','F142','F143','F144','F145','F146'};
% 'CEC2020' %#ok<NOPTS>
% FunctionSeq       = {'F1',  'F2',  'F3',  'F4',  'F5',  'F6',  'F7',  'F8',  'F9',  'F10'};
% Function_name_all = {'F147','F148','F149','F150','F151','F152','F153','F154','F155','F156'};
% 'CEC2022' %#ok<NOPTS>
% FunctionSeq       = {'F1',  'F2',  'F3',  'F4',  'F5',  'F6',  'F7',  'F8',  'F9',  'F10', 'F11', 'F12'};
% Function_name_all = {'F157','F158','F159','F160','F161','F162','F163','F164','F165','F166','F167','F168'} %#ok<NOPTS>
% 'F1-F4' %#ok<NOPTS>
% FunctionSeq       = {'F1',  'F3',  'F4'};
% Function_name_all = {'F107','F109','F110'};
% 'F5-F11' %#ok<NOPTS>
% FunctionSeq       = {'F5',  'F6',  'F7',  'F8',  'F9',  'F10', 'F11'};
% Function_name_all = {'F111','F112','F113','F114','F115','F116','F117'} %#ok<NOPTS>
% 'F12-F20' %#ok<NOPTS>
% FunctionSeq       = {'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18',  'F19',  'F20'};
% Function_name_all = {'F118','F119','F120','F121','F122','F123','F124', 'F125', 'F126'};
% 'F23-F29' %#ok<NOPTS>
% FunctionSeq       = {'F23',  'F24',  'F25',  'F26',  'F27',  'F28', 'F29'};
% Function_name_all = {'F129', 'F130', 'F131', 'F132', 'F133', 'F134','F135'};
% 'Random selected' %#ok<NOPTS>
% FunctionSeq       = {'F1',   'F3',   'F5',   'F8',   'F12',  'F17',  'F18',  'F23',  'F26',  'F29'};
% Function_name_all = {'F107', 'F109', 'F111', 'F114', 'F119', 'F124', 'F125', 'F130', 'F133', 'F136'};

%% Global parameters
SearchAgents_no=30; MaxFEs = 30 * 5000 * 2;
% algorithmName = {'MSAO', 'AO', 'CLPSO', 'MGO', 'AO', 'RIME'} %#ok<NOPTS>
% algorithmName = {'DCS2', 'DCS'} %#ok<NOPTS>
% algorithmName = {'HSA_final_0306', 'AO'} %#ok<NOPTS>
% algorithmName = {'G2', 'G22', 'G4', 'G5', 'GO2'} %#ok<NOPTS>
% algorithmName = {'SBO', 'MGO'} %#ok<NOPTS>
% algorithmName = {'SBO', 'BSA'} %#ok<NOPTS>
% algorithmName = {'QAGO', 'GO', 'AO'} %#ok<NOPTS>
% algorithmName = {'SBO', 'QAGO', 'BSA'} %#ok<NOPTS>
% algorithmName = {'SBO', 'SBO_EPM_Strategy1', 'SBO_EPM_Strategy2'} %#ok<NOPTS>
% algorithmName = {'PLO0', 'PLO_EPM_Strategy1', 'PLO_EPM_Strategy2', 'XXPLO'} %#ok<NOPTS>
% algorithmName = {'PLO0', 'XXPLO_PB_RTR'} %#ok<NOPTS>
% algorithmName = {'PLO0', 'XXPLO_PB_RTR_A', 'XXPLO_PB_RTR_ACC', 'XXPLO_PB_RTR_k', 'XXPLO_PB_RTR_FDB'} %#ok<NOPTS>
% algorithmName = {'PLO0', 'XXPLO_PB_RTR_A'} %#ok<NOPTS>
% algorithmName = {'PLO0','XXPLO', 'XXPLO_variant1_JADE', 'XXPLO_variant2_EPM'};% AvgFitnessRank:
% AvgRank:
% PLO0 : 2.4138
% XXPLO : 2.9655
% XXPLO_variant1_JADE : 1.9655
% XXPLO_variant2_EPM : 2.6552
% Rank:
% PLO0 : 2
% XXPLO : 4
% XXPLO_variant1_JADE : 1
% XXPLO_variant2_EPM : 3
% Rank:
% PLO0 : 2
% XXPLO_PB_RTR_A : 1
% XXPLO_PB_RTR_ACC : 3
% XXPLO_PB_RTR_k : 5
% XXPLO_PB_RTR_FDB : 4

% algorithmName = {'PLO0', 'XXPLO_variant1_JADE'};% AvgFitnessRank:
% AvgFitnessRank:
% PLO0 : 1  2  1  2  2  2  2  2  1  2  2  2  2  2  1  2  2  2  2  1  1  1  2  1  1  1  2  1  2
% XXPLO_variant1_JADE : 2  1  2  1  1  1  1  1  2  1  1  1  1  1  2  1  1  1  1  2  2  2  1  2  2  2  1  2  1
% AvgFitnessRankSum:
% PLO0 : 47
% XXPLO_variant1_JADE : 40
% AvgRank:
% PLO0 : 1.6207
% XXPLO_variant1_JADE : 1.3793
% Rank:
% PLO0 : 2
% XXPLO_variant1_JADE : 1

% algorithmName = {'PLO0', 'XXPLO_variant1_JADE', 'XXPLO_variant1_JADE_EPM'};% AvgFitnessRank:
% AvgRank:
% PLO0 : 2.3448
% XXPLO_variant1_JADE : 2.069
% XXPLO_variant1_JADE_EPM : 1.5862
% Rank:
% PLO0 : 3
% XXPLO_variant1_JADE : 2
% XXPLO_variant1_JADE_EPM : 1

% algorithmName = {'PLO0', 'XXPLO_variant1_EPM'};% AvgFitnessRank:
% AvgFitnessRankSum:
% PLO0 : 52
% XXPLO_variant1_EPM : 35
% AvgRank:
% PLO0 : 1.7931
% XXPLO_variant1_EPM : 1.2069
% Rank:
% PLO0 : 2
% XXPLO_variant1_EPM : 1

% algorithmName = {'XXPLO_variant1_JADE_EPM', 'XXPLO_variant1_JADE', 'XXPLO_variant1_EPM', 'PLO0'};% AvgFitnessRank:
% algorithmName = {'XXPLO_variant1_JADE_EPM', 'XXPLO_variant1_JADE'};% AvgFitnessRank:
% algorithmName = {'XXPLO_variant1_EPM', 'XXPLO_variant1_EPM_JADE_v2'};% AvgFitnessRank:
% AvgFitnessRankSum:
% XXPLO_variant1_EPM : 47
% XXPLO_variant1_EPM_JADE_v2 : 40
% AvgRank:
% XXPLO_variant1_EPM : 1.6207
% XXPLO_variant1_EPM_JADE_v2 : 1.3793
% Rank:
% XXPLO_variant1_EPM : 2
% XXPLO_variant1_EPM_JADE_v2 : 1

% algorithmName = {'XXPLO_variant1_EPM_JADE_v2', 'XXPLO_variant1_EPM', 'XXPLO_variant1_JADE_v2', 'PLO'};% AvgFitnessRank:
% AvgRank:
% XXPLO_variant1_EPM_JADE_v2 : 2.1379
% XXPLO_variant1_EPM : 2.2069
% XXPLO_variant1_JADE_v2 : 2.2759
% PLO0 : 3.3793
% Rank:
% XXPLO_variant1_EPM_JADE_v2 : 1
% XXPLO_variant1_EPM : 2
% XXPLO_variant1_JADE_v2 : 3
% PLO0 : 4

% algorithmName = {'PLOJF', 'PLOJ', 'PLOF', 'PLO'};% AvgFitnessRank:
% algorithmName = {'RIME', 'ARIME', 'ARIME2', 'ARIME3', 'ARIME4', 'ARIME5'};% AvgFitnessRank:
% algorithmName = {'RIME','ARIME4', 'GO'};% AvgFitnessRank:
% algorithmName = {'RIME','RIME_GO'};% AvgFitnessRank:
% algorithmName = {'RIME_onlySoft','RIME_onlyPunct'};% AvgFitnessRank:
% algorithmName = {'GO_onlyLearn','GO_onlyReflect', 'GO_noRandomAccept', 'GO_noMutation'};% AvgFitnessRank:
% algorithmName = {'RIME','RIME_plusGOReflect','RIME_plusGOLearnReflect','GO_noRandomAccept'};% AvgFitnessRank:
% algorithmName = {'RIME','RIME_plusGOReflect','RIME_plusGOLearnReflect','RIME_plusGOLearnReflect_SG','GO_noRandomAccept'};% AvgFitnessRank:
% algorithmName = {'RIME','RIME_plusGOReflect','RIME_plusGOLearnReflect','RIME_plusGOLearn'};% AvgFitnessRank:
algorithmName = {'LRRIME','LRIME','RRIME','RIME'};% AvgFitnessRank:

algorithmName = {'IVY','AO','DCS','YKAO'};% AvgFitnessRank:
algorithmName = {'SGAO','YKAO','YKAO_re','IVYAO'};% AvgFitnessRank:
algorithmName = {'SGAO','SAO','GAO','AO'};% AvgFitnessRank:

% algorithmName = {'MDPLO', 'CMAES'} %#ok<NOPTS>
% algorithmName = {'CGPLO_DR', 'PLO'} %#ok<NOPTS>
% algorithmName = {'CDPLO', 'CPLO', 'DPLO', 'PLO'} %#ok<NOPTS>
% algorithmName = {'CDPLO', 'CPLO'} %#ok<NOPTS>
% algorithmName = {'SBO', 'SBO0217', 'BSA'} %#ok<NOPTS>
algorithmName = {'BSA', 'BSA'} %#ok<NOPTS>
% algorithmName = {'MSAO', 'GGAO', 'TSAO', 'AO'} %#ok<NOPTS>
% algorithmName = {'MSAO', 'ALCPSO', 'BLPSO'} %#ok<NOPTS>
% algorithmName = {'MSAO', 'ALCPSO', 'BLPSO', 'CMAES', 'EOBLSSA', 'EPSDE',...
%     'GOTLBO', 'CBA', 'CLPSO', 'CAGWO', 'GWOCMA', 'RUN', 'CCMWOA', 'RDWOA', ...
%     'SMA', 'RIME', 'HGS', 'INFO', 'PLO'} %#ok<NOPTS>
% algorithmName = {'IVY', 'FLA', 'RIME', 'BKA', 'MVO', 'EOBLSSA', 'CCMWOA', 'CAGWO', ...
%     'GWOCMA', 'RDWOA', 'CBA', 'CLPSO', 'RIME'} %#ok<NOPTS>

Fold=4;
NumofRecord = 40;
foldNum=size(algorithmName,2);
funcNum=size(Function_name_all,2);
bestFitnessAll = inf * ones(foldNum, Fold, funcNum);

%% Prepare the folder
expTime = datetime();
[year, month, day] = ymd(expTime);
[hour, minute, ~] = hms(expTime);
dayStr = [num2str(year),'-',num2str(month),'-',num2str(day)];
timeStr = [num2str(hour),'_',num2str(minute)];
% 为避免出现 "WangJianFolder/WangJianFolder/expResult" 的重复目录问题，
% 将结果目录基于当前脚本所在目录构建
scriptDir = fileparts(mfilename('fullpath'));
dirName = fullfile(scriptDir, 'expResult', dayStr, [algorithmName{1} '-' timeStr]); % Store the figures
if ~exist(dirName, 'dir'); mkdir(dirName); end

%% Main process to bench
for func = 1:funcNum
    dim = 30; % For CEC17
    Function_name=Function_name_all{func};
    FunctionSeq_name=FunctionSeq{func};
    [lb,ub,dim,fhd]=Get_Functions(Function_name, dim); % For CEC17
    % [lb,ub,dim,fhd]=Get_Functions(Function_name); % For CEC22
    disp(['----------------',Function_name,'-',FunctionSeq_name,'--------------------']);
    Function_Title_name=FunctionSeq_name; % Test using
    % Function_Title_name=['F',num2str(func)]; % Experiment using
    
    %% Prepare the cell array to hold function handles
    alg_fhds = cell(1, foldNum);
    for algo = 1:foldNum
        alg_fhds{algo} = str2func(algorithmName{algo});
    end

    %% This used to store the best fitness of every run
    cg_curves=zeros(foldNum,Fold,NumofRecord);
    parfor fold=1:Fold
        display(['fold',num2str(fold)]);
        bestFitnessLocal = zeros(foldNum, 1);
        for algo=1:foldNum
            [best_pos,cg_curve]=alg_fhds{algo}(SearchAgents_no,MaxFEs,lb,ub,dim,fhd);
            if size(best_pos, 2) ~= dim
                best_pos = best_pos.Position;
            end
            bestFitness = cg_curve(end);
            bestFitnessLocal(algo, 1) = bestFitness;
            cg_curves(algo,fold,:) = MySampling(cg_curve,NumofRecord);
        end
        bestFitnessAll(:, fold, func) = bestFitnessLocal;
    end
    
    %% Plot the convergence curve
    all_curves=zeros(foldNum * Fold, NumofRecord);
    for algo = 1:foldNum
        all_curves((algo-1)*Fold+1:(algo-1)*Fold+Fold,:)=cg_curves(algo,:,:);
    end
    clf
    % figure
    screenSize = get(0, 'ScreenSize');
    figureWidth = 1000;
    figureHeight = 600;    
    % figureWidth = 500;
    % figureHeight = 300;
    xPosition = screenSize(3) - figureWidth;
    yPosition = 0;
    set(gcf,'Position',[xPosition,yPosition,figureWidth,figureHeight])
    YValue = zeros(foldNum, NumofRecord);
    for algo = 1:foldNum
        YValue(algo,:) = mean(all_curves((algo-1)*Fold+1:(algo-1)*Fold+Fold,:));
    end
    XValue = [1:NumofRecord] * (MaxFEs/NumofRecord);
    basic_linestyles = cellstr(char('-',':','-.','--'));
%     basic_Markers    = cellstr(char('o','x','+','*','s','d','v','^','<','>','p','h','.'));
    basic_Markers    = cellstr(char(''));
    nLines = foldNum;
    MarkerEdgeColors = hsv(nLines);
    linestyles = repmat(basic_linestyles,ceil(nLines/numel(basic_linestyles)),1);
    Markers = repmat(basic_Markers,ceil(nLines/numel(basic_Markers)),1);
    set(gcf,'color','white')
    for algo=1:foldNum
        if algo == 1
            LineWidthSet = 1.5;
        else
            LineWidthSet = 1.5;
        end
        semilogy(XValue,YValue(algo,:),[linestyles{algo} Markers{algo}],'LineWidth', LineWidthSet, ...
            'Color',MarkerEdgeColors(algo,:), 'DisplayName', algorithmName{algo});
        hold on;
    end
    hold off;
    xlabel('FEs');
    ylabel('Best fitness');
    title(Function_Title_name);
    % legend show;
    legend('show', 'Box', 'off');
    % legend boxoff;
    %% Beautify the figure
    a=findobj(gcf); 
    allaxes=findall(a,'Type','axes');
    alltext=findall(a,'Type','text');
    set(allaxes,'FontName','Times','LineWidth',1,'FontSize',13,'FontWeight','bold');
    set(alltext,'FontName','Times','FontSize',13,'FontWeight','bold')
    set(gcf, 'PaperUnits', 'inches');
    krare=2.5;
    figureWidth=krare*5/3 ;
    figureHeight=krare*4/3;
    set(gcf, 'PaperPosition', [0,0,figureWidth,figureHeight]); 
    set(gca, ... % Beautify axes
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );
    axis tight;
    CCFigName = [dirName, '/', algorithmName{1}, '-', Function_Title_name, '-CC'];
    print(CCFigName,'-dtiff', '-r300')
    saveas(gcf, CCFigName, 'fig')
end
close

%% Take the fold average and compress to two-dimensional array
AvgFitness = mean(bestFitnessAll, 2);
AvgFitness = squeeze(AvgFitness);
% disp(AvgFitness);

%% Take the fold average and compress to two-dimensional array
StdFitness = std(bestFitnessAll, 0, 2);
StdFitness = squeeze(StdFitness);
% disp(StdFitness);

%% Rank every column (based on same function sort)
AvgFitnessRank = AvgFitness;
for col = 1:size(AvgFitness, 2)
    [~, ~, ColRank] = unique(AvgFitness(:, col));
    AvgFitnessRank(:, col) = ColRank;
end
AvgFitnessRankSum = sum(AvgFitnessRank, 2);
[~, ~, AlgoRank] = unique(AvgFitnessRankSum);

%% Show the rank of algorithms in each functions
disp('AvgFitnessRank:');
for algo = 1:foldNum
    disp([algorithmName{algo}, ' : ', num2str(AvgFitnessRank(algo, :))]);
end

%% Show the rank sum of all test functions
disp('AvgFitnessRankSum:');
for algo = 1:foldNum
    disp([algorithmName{algo}, ' : ', num2str(AvgFitnessRankSum(algo))]);
end

%% Show the average rank
disp('AvgRank:');
for algo = 1:foldNum
    disp([algorithmName{algo}, ' : ', num2str(AvgFitnessRankSum(algo) ./ funcNum)]);
end

%% Show the integer rank
Rank = zeros(1, foldNum);
[~, index] = sort(AvgFitnessRankSum ./ funcNum);
for algo = 1:algo
    Rank(index(algo)) = algo;
end
disp('Rank:');
for algo = 1:foldNum
    disp([algorithmName{algo}, ' : ', num2str(Rank(algo))]);
end

%% To write the avg and std data to excel file
fileName = [dirName, '/', algorithmName{1}];
xlsFileName = [fileName, '.xlsx'];
% To prepare the result data and write to xls file
result = {'Func', 'Metric'};
result{2*funcNum+2, 1} = 'ARV';
result{2*funcNum+3, 1} = 'Rank';
for algo = 1:foldNum
    result{1, algo+2} = algorithmName{algo};
    result{2*funcNum+2, algo+2} = num2str(AvgFitnessRankSum(algo) ./ funcNum);
    result{2*funcNum+3, algo+2} = num2str(Rank(algo));
    for func = 1:funcNum
%         if isempty(result{2*func, 1}), result{2*func, 1} = ['F', num2str(func)]; end
%         if isempty(result{2*func, 2}), result{2*func, 2} = 'Avg'; end
%         if isempty(result{2*func+1, 2}), result{2*func+1, 2} = 'Std'; end
        result{2*func, 1} = ['F', num2str(func)];
        result{2*func, 2} = 'Avg';
        result{2*func+1, 2} = 'Std';
        result{2*func, algo+2} = AvgFitness(algo, func);
        result{2*func+1, algo+2} = StdFitness(algo, func);
    end
end
writecell(result, xlsFileName, 'sheet', 'cmpResult');

%% To write the pValue data to excel file
value = {'Func'};
for algo = 1:foldNum-1
    value{1, algo+1} = algorithmName{algo+1}; % Write the first row
end
for func = 1:funcNum
    value{func+1, 1} = ['F', num2str(func)]; % Write the first column
end
% Now write the pValue data
for func = 1:funcNum
    for algo = 1:foldNum - 1
        value{func+1, algo+1} = num2str(signrank(bestFitnessAll(1, :, func), bestFitnessAll(algo+1, :, func)));
    end
end
writecell(value, xlsFileName, 'sheet', 'pValue');

%% To write the rank & pValue excel file
RankPValue{1, 1} = 'Func';
RankPValue{1, 2} = algorithmName{1};
for func = 1:funcNum % The first column
    RankPValue{func+1, 1} = ['F', num2str(func)];
end
for algo = 1:foldNum-1 % The first row
    RankPValue{1, 2*algo+1} = algorithmName{algo+1};
    RankPValue{1, 2*(algo+1)} = 'pValue';
end
for func = 1:funcNum % The algorithmName{1} rank
    RankPValue{func+1, 2} = num2str(AvgFitnessRank(1, func));
end
for algo = 1:foldNum-1 % The data table
    for func = 1:funcNum % The other algorithm's rank
        RankPValue{func+1, 2*algo+1} = num2str(AvgFitnessRank(algo+1, func));
        RankPValue{func+1, 2*algo+2} = value{func+1, algo+1};
    end
end
writecell(RankPValue, xlsFileName, 'sheet', 'rank & pValue');

%% Write the pValue & sign table
pValueSign{1, 1} = 'Func';
for func = 1:funcNum % The first column
    pValueSign{func+1, 1} = ['F', num2str(func)];
end
for algo = 1:foldNum-1 % The first row
    pValueSign{1, algo+1} = algorithmName{algo+1};
end
% RankPValue = xlsread(xlsFileName, 'rank & pValue');
RankPValue = readmatrix(xlsFileName, 'sheet', 'rank & pValue');
RankPValue = RankPValue(:, 2:end);
pValueSign = cell(funcNum+2, foldNum);
pValueSign{1, 1} = 'Func';
pValueSign{funcNum+2, 1} = '+/-/=';
pValueSign{funcNum+3, 1} = '+/-/=';
for func = 1:funcNum % The first column
    pValueSign{func+1, 1} = ['F', num2str(func)];
end
for algo = 1:foldNum-1 % The first row
    pValueSign{1, 2*algo} = algorithmName{algo+1};
end
for algo = 1:foldNum-1 % The pValue and sign table
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
        pValueSign{func+1, 2*algo} = value{func+1, algo+1};
    end
    cntEqual = funcNum - cntBetter - cntWorse;
    pValueSign{func+2, 2*algo} = [num2str(cntBetter), '/', num2str(cntWorse), ...
        '/', num2str(cntEqual)];
    pValueSign{func+3, algo+1} = [num2str(cntBetter), '/', num2str(cntWorse), ...
        '/', num2str(cntEqual)];
end
writecell(pValueSign, xlsFileName, 'sheet', 'pValue & sign');

%% To write the fold final fitness for boxplot
for func = 1:funcNum
    data = cell(Fold, foldNum);
    data{1, 1} = 'Fold';
    for algo = 1:foldNum
        data{1, algo+1} = algorithmName{algo};
    end
    for algo = 1:foldNum
        for fold = 1:Fold
            data{fold+1, 1} = ['fold', num2str(fold)];
            data{fold+1, algo+1} = bestFitnessAll(algo, fold, func);
        end
    end
    sheetName = ['F', num2str(func)];
    writecell(data, xlsFileName, 'sheet', sheetName);
end

% %% Show the rank of the algorithms of all test functions
% disp('AlgoRank:');
% disp(AlgoRank);

% %% Plot every algorithm's rank in each function as a line
% figure;
% hold on;
% for i = 1:size(AvgFitnessRank, 1)
%     plot(1:size(AvgFitnessRank, 2), AvgFitnessRank(i, :), 'DisplayName', algorithmName{i});
% end
% hold off;
% % Add label and legend
% xticks(1:size(AvgFitnessRank, 2));
% xlabel('Function Index');
% ylabel('Rank');
% legend show;
% title('Plot of every function rank over functions');、

%% To Plot the boxplot
for func = 1:funcNum
    figureTitle = ['F', num2str(func)];
    plotData = readmatrix(xlsFileName, 'sheet', figureTitle);
    plotData = plotData(:, 2:end);
    fenceMargin = 0.14;
    rectangleMargin = 0.28;
    outMarker = '.';
    outMarkerSize = 3;
    outMarkerEdgeColor = [0.6 0.6 0.6];
    outMarkerFaceColor = [rand rand rand];
    algoNum = size(plotData, 2);
    
    colors = slanCM(100, algoNum);

    xlim([0.5 algoNum+0.5]);

    for algo=1:algoNum
    
        % Calculate the mean and intervals
        boxData = sort(plotData(:, algo));
        quater2 = median(boxData);
        quater1 = median(boxData(boxData <= quater2));
        quater3 = median(boxData(boxData >= quater2));
        IQR = quater3 - quater1;
        fenceLow = quater1 - 1.5 * IQR;
        fenceUp = quater3 + 1.5 * IQR;
        % fenceLow = min(boxData(boxData >= quater1-1.5*IQR));
        % fenceUp = max(boxData(boxData <= quater3+1.5*IQR));
        % outlierLow = boxData(boxData < quater1-1.5*IQR);
        % outlierUp = boxData(boxData > quater3+1.5*IQR);
        meanDot = mean(plotData(:,algo));

        % hold on;
        % plot((algo).*ones(size(outlierUp)),outlierUp,...
        %     'LineStyle','none',...
        %     'Marker',outMarker,...
        %     'MarkerSize',outMarkerSize,...
        %     'MarkerEdgeColor',outMarkerEdgeColor,...
        %     'MarkerFaceColor',outMarkerFaceColor);
        % plot((algo).*ones(size(outlierLow)),outlierLow,...
        %     'LineStyle','none',...
        %     'Marker',outMarker,...
        %     'MarkerSize',outMarkerSize,...
        %     'MarkerEdgeColor',outMarkerEdgeColor,...
        %     'MarkerFaceColor',outMarkerFaceColor);

        if sum(isnan(boxData)) ~= 0
            disp('DEBUG: NaN number exist in boxData');
        end
       
        if meanDot < fenceLow || meanDot > fenceUp
            disp('DEBUG: meanDot out of range of fenceUp and fenceLow');
        end
        
        hold on;
        % Plot fence
        line([algo algo],[fenceUp fenceLow],...
            'Color','k','LineStyle',':');
        line([algo-fenceMargin algo+fenceMargin],[fenceUp fenceUp],...
            'Color','k');
        line([algo-fenceMargin algo+fenceMargin],[fenceLow fenceLow],...
            'Color','k');
        
        % Plot quantile
        if quater3 > quater1
            rectangle('Position',[algo-rectangleMargin quater1 2*rectangleMargin quater3-quater1],...
                'EdgeColor','k','FaceColor',colors(algo, :));
        end
        
        % Plot median
        line([algo-rectangleMargin algo+rectangleMargin],[quater2 quater2],...
            'Color','k','LineWidth',1);
        
        % Plot datapoint
        hold on;
        % plot(algo.*ones(size(plotData, 1)), plotData(:, algo)', ...
        % Plot meanDot
        plot(algo, meanDot, ...
            'LineStyle','none', ...
            'Marker','o', ...
            'MarkerEdgeColor','k', ...
            'MarkerFaceColor',colors(algo, :));
        hold off;
    end
    
    box on;
    
    set(gcf, 'PaperUnits', 'inches');
    set(gca,'XTick',1:algoNum);
    set(gca,'XTickLabel',algorithmName);
    %% Image beautification
    a=findobj(gcf);
    allaxes=findall(a,'Type','axes');
    alltext=findall(a,'Type','text');
    set(allaxes,'FontName','Times','LineWidth',1,'FontSize',12,'FontWeight','bold');
    set(alltext,'FontName','Times','FontSize',12,'FontWeight','bold')
    set(gcf, 'PaperUnits', 'inches');
    krare=2.5;
    figureWidth=krare*5/3 ;
    figureHeight=krare*4/3;
    set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);
    set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XGrid'       , 'off'      , ...
    'GridColor'   , [1.0 1.0 1.0],...    
    'MinorGridColor',[1.0 1.0 1.0],...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );
    % axis tight
    grid off
    box on
    % xlabel('Algorithms'); % Set the X-axis label
    ylabel('Fitness');
    title(figureTitle);
    boxFigureName = [dirName, '/', algorithmName{1}, '-', figureTitle, '-boxplot'];
    print(boxFigureName, '-dtiff', '-r300');
    saveas(gcf, boxFigureName, 'fig')

    ylim('auto')
    fontSize = 10;
    labelX = findobj(gca, 'Type', 'text');
    rotation = 90;
    for cnt = 1:length(labelX)
        set(labelX(cnt),    'FontSize', fontSize,...
            'Rotation', rotation, ...
            'String', tickLabelStr{algo}, ...
            'HorizontalAlignment', 'right',...
            'FontWeight','bold',...
            'interpreter','latex');
        set(gcf, 'PaperUnits', 'inches');
    end
    close all
end
toc
% end