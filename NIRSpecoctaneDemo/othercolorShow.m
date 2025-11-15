function othercolorShow(colors,id0,num)
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);

% Background is white
set(gcf,'color','white');  

for k = id0:num+id0-1
       subplot(ceil(num/10),10, k-id0+1);
       c = othercolor(k);
       imagesc(reshape(c,1,size(c,1),size(c,2)));
       title(['No.',num2str(k),': ',char(colors(k))],'interpreter','none');
       axis off;
end
end