clear all 
clc
% add function
addpath([pwd '/CCMWOA'])
addpath([pwd '/Origin'])
addpath([pwd '/subFunction'])
addpath([pwd '/CCMWOA/ALCPSO'])
SearchAgents_no=20; % Number of search agents
%1,2,3,4,5,6,7,8,10,11
Function_name='F24'; % Name of the test function that can be from F1 to F23 (F8,)
Max_iteration=100; % Maximum numbef of iterations 
dim=10;
MaxFEs=dim*10000;
for funcNum=3:4
    Function_name=['F',num2str(funcNum)];
    % Load details of the selected benchmark function
    [lb,ub,dim,fobj]=Get_Functions(Function_name,10);
    [~,cg_curve1]=BA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj);
end

semilogy(cg_curve1,'Color',[1,0,0],'LineWidth', 1.5);
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on


        



