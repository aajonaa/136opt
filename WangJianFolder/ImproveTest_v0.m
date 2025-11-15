%% Initial setting
clear; close all;  clc; tic;
addpath(genpath(pwd));

%% Test functions choose
'CEC2017' %#ok<NOPTS>
FunctionSeq       = {'F1',  'F3',  'F4',  'F5',  'F6',  'F7',  'F8',  'F9',  'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'F17', 'F18', 'F19', 'F20', 'F21', 'F22', 'F23', 'F24', 'F25', 'F26', 'F27', 'F28', 'F29', 'F30'};
Function_name_all = {'F107','F109','F110','F111','F112','F113','F114','F115','F116','F117','F118','F119','F120','F121','F122','F123','F124','F125','F126','F127','F128','F129','F130','F131','F132','F133','F134','F135','F136'};
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
algorithmName = {'AO', 'AO0', 'CLPSO'} %#ok<NOPTS>
Flod=4;
NumofRecord = 40;
algoNum=size(algorithmName,2);
funcNum=size(Function_name_all,2);
bestFitnessAll = inf * ones(algoNum, Flod, funcNum);

%% Prepare the folder
daystr = datestr(now, 'mm-dd');
timestr = datestr(now, 'mm-dd_HH_MM');
dirName = ['WangJianFolder/expResult/', daystr]; % Store the figures
dirName = [dirName, '/', algorithmName{1}, '-', timestr]; % Store the figures
mkdir(dirName);

%% Main process to bench
for func = 1:funcNum
    dim = 30;
    Function_name=Function_name_all{func};
    FunctionSeq_name=FunctionSeq{func};
    [lb,ub,dim,fhd]=Get_Functions(Function_name, dim);
    disp(['----------------',Function_name,'-',FunctionSeq_name,'--------------------']);
%     Function_Title_name=FunctionSeq_name; % Test using
    Function_Title_name=['F',num2str(func)]; % Experiment using
    
    %% Prepare the cell array to hold function handles
    alg_fhds = cell(1, algoNum);
    for algo = 1:algoNum
        alg_fhds{algo} = str2func(algorithmName{algo});
    end

    %% This used to store the best fitness of every run
    cg_curves=zeros(algoNum,Flod,NumofRecord);
    parfor flod=1:Flod
        display(['flod',num2str(flod)]);
        bestFitnessLocal = zeros(algoNum, 1);
        for algo=1:algoNum
            [best_pos,cg_curve]=alg_fhds{algo}(SearchAgents_no,MaxFEs,lb,ub,dim,fhd);
            bestFitness = fhd(best_pos);
            bestFitnessLocal(algo, 1) = bestFitness;
            cg_curves(algo,flod,:) = MySampling(cg_curve,NumofRecord);
        end
        bestFitnessAll(:, flod, func) = bestFitnessLocal;
    end
    
    %% Plot the convergence curve
    all_curves=zeros(algoNum * Flod, NumofRecord);
    for algo = 1:algoNum
        all_curves((algo-1)*Flod+1:(algo-1)*Flod+Flod,:)=cg_curves(algo,:,:);
    end
%     clf
    figure
    screenSize = get(0, 'ScreenSize');
    figureWidth = 1000;
    figureHeight = 600;
    xPosition = screenSize(3) - figureWidth;
    yPosition = 0;
    set(gcf,'Position',[xPosition,yPosition,figureWidth,figureHeight])
    YValue = zeros(algoNum, NumofRecord);
    for algo = 1:algoNum
        YValue(algo,:) = mean(all_curves((algo-1)*Flod+1:(algo-1)*Flod+Flod,:));
    end
    XValue = [1:NumofRecord] * (MaxFEs/NumofRecord);
    basic_linestyles = cellstr(char('-',':','-.','--'));
%     basic_Markers    = cellstr(char('o','x','+','*','s','d','v','^','<','>','p','h','.'));
    basic_Markers    = cellstr(char(''));
    nLines = algoNum;
    MarkerEdgeColors = hsv(nLines);
    linestyles = repmat(basic_linestyles,ceil(nLines/numel(basic_linestyles)),1);
    Markers = repmat(basic_Markers,ceil(nLines/numel(basic_Markers)),1);
    set(gcf,'color','white')
    for algo=1:algoNum
        if algo == 1
            LineWidthSet = 3;
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
%     legend show;
%     legend boxoff;
    % Set the figure's property (Beautify)
    a=findobj(gcf); 
    allaxes=findall(a,'Type','axes');
    alltext=findall(a,'Type','text');
    set(allaxes,'FontName','Times','LineWidth',1,'FontSize',13,'FontWeight','bold');
    set(alltext,'FontName','Times','FontSize',13,'FontWeight','bold')
    set(gcf, 'PaperUnits', 'inches');
    krare=2.5;
    figureWidth=krare*5/3 ;
    figureHeight=krare*4/3;
    set(gcf, 'PaperPosition', [xPosition,yPosition,figureWidth,figureHeight]); 
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
    figName = [dirName, '/', algorithmName{1}, '-', Function_Title_name];
    print(figName,'-dtiff', '-r300')
    saveas(gcf, figName, 'fig')
end

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
for algo = 1:algoNum
    disp([algorithmName{algo}, ' : ', num2str(AvgFitnessRank(algo, :))]);
end

%% Show the rank sum of all test functions
disp('AvgFitnessRankSum:');
for algo = 1:algoNum
    disp([algorithmName{algo}, ' : ', num2str(AvgFitnessRankSum(algo))]);
end

%% Show the average rank
disp('AvgRank:');
for algo = 1:algoNum
    disp([algorithmName{algo}, ' : ', num2str(AvgFitnessRankSum(algo) ./ funcNum)]);
end

%% Show the integer rank
Rank = zeros(1, algoNum);
[~, index] = sort(AvgFitnessRankSum ./ funcNum);
for i = 1:algo
    Rank(index(i)) = i;
end
disp('Rank:');
for algo = 1:algoNum
    disp([algorithmName{algo}, ' : ', num2str(Rank(algo))]);
end

%% To write the avg and std data to excel file
fileName = [dirName, '/', algorithmName{1}, '-', timestr];
xlsFileName = [fileName, '.xlsx'];
% To prepare the result data and write to xls file
result = {'Func', 'Metric'};
result{2*funcNum+2, 1} = 'ARV';
result{2*funcNum+3, 1} = 'Rank';
for algo = 1:algoNum
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
xlswrite(xlsFileName, result, 'cmpResult');

%% To write the pValue data to excel file
value = {'Func'};
for algo = 1:algoNum-1
    value{1, algo+1} = algorithmName{algo+1}; % Write the first row
end
for func = 1:funcNum
    value{func+1, 1} = ['F', num2str(func)]; % Write the first column
end
% Now write the pValue data
for func = 1:funcNum
    for algo = 1:algoNum - 1
        value{func+1, algo+1} = signrank(bestFitnessAll(1, :, func), bestFitnessAll(algo+1, :, func));
    end
end
xlswrite(xlsFileName, value, 'pValue');

%% To write the rank & pValue excel file
RankPValue{1, 1} = 'Func';
RankPValue{1, 2} = algorithmName{1};
for func = 1:funcNum % The first column
    RankPValue{func+1, 1} = ['F', num2str(func)];
end
for algo = 1:algoNum-1 % The first row
    RankPValue{1, 2*algo+1} = algorithmName{algo+1};
    RankPValue{1, 2*(algo+1)} = 'pValue';
end
for func = 1:funcNum % The algorithmName{1} rank
    RankPValue{func+1, 2} = num2str(AvgFitnessRank(1, func));
end
for algo = 1:algoNum-1 % The data table
    for func = 1:funcNum % The other algorithm's rank
        RankPValue{func+1, 2*algo+1} = num2str(AvgFitnessRank(algo+1, func));
        RankPValue{func+1, 2*algo+2} = value{func+1, algo+1};
    end
end
xlswrite(xlsFileName, RankPValue, 'rank & pValue');

%% Write the pValue & sign table
pValueSign{1, 1} = 'Func';
for func = 1:funcNum % The first column
    pValueSign{func+1, 1} = ['F', num2str(func)];
end
for algo = 1:algoNum-1 % The first row
    pValueSign{1, algo+1} = algorithmName{algo+1};
end
RankPValue = xlsread(xlsFileName, 'rank & pValue');
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
        pValueSign{func+1, 2*algo} = value{func+1, algo+1};
    end
    cntEqual = funcNum - cntBetter - cntWorse;
    pValueSign{func+2, 2*algo} = [num2str(cntBetter), '/', num2str(cntWorse), ...
        '/', num2str(cntEqual)];
    pValueSign{func+3, algo+1} = [num2str(cntBetter), '/', num2str(cntWorse), ...
        '/', num2str(cntEqual)];
end
xlswrite(xlsFileName, pValueSign, 'pValue & sign');

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
% title('Plot of every function rank over functions');

toc;
% close all;