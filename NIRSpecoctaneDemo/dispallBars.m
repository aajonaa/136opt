function dispallBars(Color,ids,ide,cmlib) 
% Display the all colorbars
close all
fig = figure('Units','normalized','OuterPosition',[0.25, 0, 0.5,1]); 
set(gcf,'Color','White');  % Background is white
y = linspace(0.94, 0.06, ide-ids+1);
i = 1;
for k = Color(ids:ide)    
    ax = axes(fig,'Position',[0.37 y(i) 0.4 0.0096]);
    switch cmlib
        case {'dip'}            
            map = idipCM(k);   %  k is a string of colorType arrays 
        case {'slan'}            
            map = slanCM(k);  %  Same as above
        case {'ncl'}
            map = nclCM(k);    %  Same as above
        otherwise
            map=[];
    end    
    imagesc(ax,(1:256));
    colormap(ax,map); 
    axis(ax,'off');
    % axis off
    text(ax,-100,1,['[',num2str(i+ids-1),']',' ',char((Color(i+ids-1)))],'Color','b',...
        'Fontname','Times New Roman',...
        'FontSize',10,'FontWeight','bold','Interpreter','none')
    i = i+1;
    pause(0.001);
end
end