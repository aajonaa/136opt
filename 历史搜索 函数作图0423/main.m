
clear all 
clc
addpath(genpath(pwd));
SearchAgents_no=100; % Number of search agents

Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Function_name='F3'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Function_name='F4'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Function_name='F6'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Function_name='F7'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Function_name='F9'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Function_name='F11'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Function_name='F15'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Function_name='F18'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

Max_iteration=500; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions(Function_name);

% [Target_score,Target_pos,GOA_cg_curve, Trajectories,fitness_history, position_history]=GOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

figure('Position',[454   445   894   297])
%Draw search space
% subplot(1,5,1);
func_plot(Function_name);
% title('Parameter space')
% xlabel('x_1');
% ylabel('x_2');
% zlabel([Function_name,'( x_1 , x_2 )'])
box on
axis tight
% meihua(gcf,Function_name)
% subplot(1,5,2);
% hold on
% for k1 = 1: size(position_history,1)
%     for k2 = 1: size(position_history,2)
%         plot(position_history(k1,k2,1),position_history(k1,k2,2),'.','markersize',1,'MarkerEdgeColor','k','markerfacecolor','k');
%     end
% end
% plot(Target_pos(1),Target_pos(2),'.','markersize',10,'MarkerEdgeColor','r','markerfacecolor','r');
% title('Search history (x1 and x2 only)')
% xlabel('x1')
% ylabel('x2')
% box on
% axis tight
% 
% subplot(1,5,3);
% hold on
% plot(Trajectories(1,:));
% title('Trajectory of 1st grasshopper')
% xlabel('Iteration#')
% box on
% axis tight
% 
% subplot(1,5,4);
% hold on
% plot(mean(fitness_history));
% title('Average fitness of all grasshoppers')
% xlabel('Iteration#')
% box on
% axis tight
% 
% %Draw objective space
% subplot(1,5,5);
% semilogy(GOA_cg_curve,'Color','r')
% title('Convergence curve')
% xlabel('Iteration#');
% ylabel('Best score obtained so far');
% box on
% axis tight
% set(gcf, 'position' , [39         479        1727         267]);
% 
% 
% display(['The best solution obtained by GOA is : ', num2str(Target_pos)]);
% display(['The best optimal value of the objective funciton found by GOA is : ', num2str(Target_score)]);
% close all

filename = 'SBO';
saveas(gcf,[filename,'-',Function_name,'-3D'],'fig')
% print(filename ,'-dtiff', '-r300'); %<-Save as PNG with 300 DPI


a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Times','LineWidth',1,'FontSize',13,'FontWeight','bold');
% set(alllines,'Linewidth',1);
set(alltext,'FontName','Times','FontSize',13,'FontWeight','bold')


% 
set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperUnits', 'centimeters');


krare=2.5;
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

% filename = 'SBO-F15-hs';
% filename = 'SBO-F18-hs';


% filename = 'SBO-F1-toa';
% filename = 'SBO-F3-toa';
% filename = 'SBO-F4-toa';
% filename = 'SBO-F6-toa';
% filename = 'SBO-F9-toa';
% filename = 'SBO-F11-toa';
% filename = 'SBO-F15-toa';
% filename = 'SBO-F18-toa';

% filename = 'SBO-F1-avgf';
% filename = 'SBO-F3-avgf';
% filename = 'SBO-F4-avgf';
% filename = 'SBO-F6-avgf';
% filename = 'SBO-F9-avgf';
% filename = 'SBO-F11-avgf';
% filename = 'SBO-F15-avgf';
% filename = 'SBO-F18-avgf';




print([filename,'-',Function_name,'-3D'], '-dtiff', '-r300'); %<-Save as PNG with 300 DPI


