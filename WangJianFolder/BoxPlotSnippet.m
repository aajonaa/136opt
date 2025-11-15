clear
clc
xlsFileName = 'D:/Jona/136个函数评估框架1_jona/WangJianFolder/expResult/2024-6-22/MSAO-4_35/MSAO.xlsx';
algorithmName = {'MSAO', 'AO', 'RIME', 'DE', 'MVO', 'EOBLSSA', 'CCMWOA', 'GOTLBO', ...
    'GWOCMA', 'RDWOA', 'CBA', 'CLPSO'} %#ok<NOPTS>
funcNum = 29;
algoNum = 12;
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
    expTime = datetime();
    [year, month, day] = ymd(expTime);
    [hour, minute, ~] = hms(expTime);
    dayStr = [num2str(year),'-',num2str(month),'-',num2str(day)];
    timeStr = [num2str(hour),'_',num2str(minute)];
    dirName = ['WangJianFolder/BoxPlot/', dayStr, ...
    '/', algorithmName{1}, '-', timeStr]; % Store the figures
    mkdir(dirName);
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