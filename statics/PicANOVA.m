function  PicANOVA(xlsfilename)
%--------------------------------------------------------------------------
%                     调用anova1函数作单因素一元方差分析
%--------------------------------------------------------------------------
clc
close all
clear
%%%%%%%%%%%%
% xlsfilename = 'OBLFA_GWO-03-05_21_22.xlsx';
readFilename=xlsfilename;     
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
%%%%%%%%%
for i=1:functioon_num
    Function_name=['F',num2str(i)]; 
%     Function_name={Function_name};
%     xlswrite(writeFilename, Function_name, sheet_name, ['A',num2str(i+1)]);   %pValue表中第一列
% end

%********************************读取数据***********************************
% 读取文件的第1个工作表中的数据
% [x,y]=xlsread('OBLFA_GWO-03-05_21_22.xlsx','F13');
[x,y]=xlsread(xlsfilename,Function_name);
% 提取矩阵x的第n列数据
score = x(:,end);
% 提取元胞数组y的第4列的第2行至最后一行数据，
method = y(1:end,1);
%*******************************方差分析************************************
[p,table,stats] = anova1(score,method)  % 单因素一元方差分析


%%%%% 美化代码

a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Times','LineWidth',1,'FontSize',12,'FontWeight','bold');
% set(alllines,'Linewidth',1);
set(alltext,'FontName','Times','FontSize',12,'FontWeight','bold')


% 
set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperUnits', 'centimeters');


krare=2;
% x_width=krare*1.618 ;
x_width=krare*5/3 ;

%  x_width=3*1;

y_width=krare*4/3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 4])

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

% axis tight
% grid on
% box on

print('-dtiff', '-r300'); %<-Save as PNG with 300 DPI

end