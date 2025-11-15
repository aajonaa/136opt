clear all %#ok<CLALL>
close all
clc

N=30; % Number of search agents

%% FATA
Function_name='F1'; % Name of the test function, range from F1-F13

FEs=300000; % Maximum number of evaluation times

dimSize = 30;   %dimension size

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_FATA(Function_name,dimSize);

[bestPositions,Convergence_curve]=FATA_FEs(N,FEs,lb,ub,dim,fobj);
Convergence_curve=MySampling(Convergence_curve,FEs);

%Draw objective space
figure,
hold on
semilogy(Convergence_curve,'Color','R','LineWidth',4);
title('Convergence curve')
xlabel('FEs');
ylabel('Best fitness obtained so far');
axis tight
grid off
box on
legend('FATA\_FEs')

display(['The best location of FATA_FEs is: ', num2str(bestPositions)]);


        



