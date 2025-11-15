str = {'hs', 'toa', 'avgf', 'carve'};
for i = 1:4
    index = 21;
    openFilename = ['HDTAO-F21-', str{i}, '.fig'];
    % openfig(openFilename);

    figHandle = openfig(openFilename, 'reuse');
    a=findobj(gcf); % get the handles associated with the current figure

    allaxes=findall(a,'Type','axes');
    set(allaxes,'FontName','Times','LineWidth',1,'FontSize',13,'FontWeight','bold');

    alltext=findall(a,'Type','text');
    set(alltext,'FontName','Times','FontSize',13,'FontWeight','bold')

    set(gcf, 'PaperUnits', 'inches');

    krare=2.5;
    x_width=krare*5/3 ;

    y_width=krare*4/3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

    set(gca, ...
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
%     filename = 'Beautified-Zoomed-Modified-RMRIME-F28-curve';
%     print(filename, '-dtiff', '-r300'); %<-Save as PNG with 300 DPI
    % saveFilename = ['Beautified-HDTAO-F', num2str(index), '-diversity'];
    saveFilename = ['Beautified-HDTAO-F', num2str(index), str{i}];
    print(saveFilename, '-dtiff', '-r300'); %<-Save as PNG with 300 DPI

    close(figHandle);
end


% %% CCBeautify_v0
% figHandle = openfig('Modified-RMRIME-F28-curve.fig', 'reuse');
% 
% a=findobj(gcf); % get the handles associated with the current figure
% 
% allaxes=findall(a,'Type','axes');
% set(allaxes,'FontName','Times','LineWidth',1,'FontSize',13,'FontWeight','bold');
% 
% alltext=findall(a,'Type','text');
% set(alltext,'FontName','Times','FontSize',13,'FontWeight','bold')
% 
% set(gcf, 'PaperUnits', 'inches');
% 
% krare=2.5;
% x_width=krare*5/3 ;
% 
% y_width=krare*4/3;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% 
% set(gca, ...
%     'Box'         , 'on'     , ...
%     'TickDir'     , 'in'     , ...
%     'TickLength'  , [.02 .02] , ...
%     'XMinorTick'  , 'on'      , ...
%     'YMinorTick'  , 'on'      , ...
%     'YGrid'       , 'off'      , ...
%     'XGrid'       , 'off'      , ...
%     'XColor'      , [.3 .3 .3], ...
%     'YColor'      , [.3 .3 .3], ...
%     'LineWidth'   , 1         );
% filename = 'Beautified-Zoomed-Modified-RMRIME-F28-curve';
% print(filename, '-dtiff', '-r300'); %<-Save as PNG with 300 DPI
% 
% close(figHandle);
