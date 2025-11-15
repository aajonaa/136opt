clear all 
clc
% add function
addpath([pwd '/subFunction'])
addpath([pwd '/Origin'])
addpath([pwd '/CCMWOA'])

SearchAgents_no=20; % Number of search agents
%1,2,3,4,5,6,7,8,10,11
Function_name='F1'; % Name of the test function 
Max_iteration=100; % Maximum numbef of iterations
MaxFEs=SearchAgents_no*1000;
% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name,10);

[~,cg_curve1]=CCMWOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj);
[~,cg_curve2]=DECLS(SearchAgents_no,MaxFEs,lb,ub,dim,fobj);

lengedstr={'DE','DECLS'};
display(['------------------------------------']);

%% Draw search space
 set(gcf, 'position', [10 10 800 400]);
% set(gcf, 'position', [10 10 300 300]);
func_plot(Function_name);
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel(['F1','( x_1 , x_2 )'])

%% Draw objective space
subplot(1,2,2);
cg_curve1=MySampling(cg_curve1,20);
semilogy(cg_curve1,'Color',[1,0,0],'LineWidth', 1.5);
hold on
cg_curve2=MySampling(cg_curve2,20);
semilogy(cg_curve2,'Color',[0,1,0],'LineWidth', 1.5)
hold on 

title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend(lengedstr);



        



