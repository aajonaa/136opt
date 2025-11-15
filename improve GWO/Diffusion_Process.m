
%This function is used to mimic diffusion process, and creates some
%new Positions based on Gaussian Walks.
%**************************************************************************
%The input function is:                                                   %
%Position: the input Position which is going to be diffused                     %
%S: structure of problem information                                      %
%g: generation number                                                     %
%BestPosition: the best Position in group                                       %                               %
%==========================================================================
%The output function is:                                                  %
%createPosition: the new Positions created by Diffusion process                 %
%fitness: the value of fitness function                                   %
%**************************************************************************
function [createPosition, fitness] = Diffusion_Process(Position,Maximum_Diffusion,Walk,lb,ub,g,BestPosition,fobj)
    %calculating the maximum diffusion for each Position
    NumDiffiusion = Maximum_Diffusion;
    New_Position = Position;
    
    %Diffiusing Part*******************************************************
    for i = 1 : NumDiffiusion
        %consider which walks should be selected.
        if rand < Walk 
            GeneratePosition = normrnd(BestPosition, (log(g)/g)*(abs((Position - BestPosition))), [1 size(Position,2)]) + ...
                (randn*BestPosition - randn*Position);
        else
            GeneratePosition = normrnd(Position, (log(g)/g)*(abs((Position - BestPosition))),...
                [1 size(Position,2)]);
        end
        New_Position = [New_Position;GeneratePosition];
    end
    %check bounds of New Position
    New_Position = Bound_Checking(New_Position,lb,ub);
    %sorting fitness
    fitness = [];
    for i = 1 : size(New_Position,1)
        g = g + 1;
        fitness = [fitness;fobj(New_Position(i,:))];
    end
    [fit_value,fit_index] = sort(fitness);
    fitness = fit_value(1,1);
    New_Position = New_Position(fit_index,:);
    createPosition = New_Position(1,:);
    %======================================================================
end

function fit = DeJong(x)
    fit=sum((x).^2);
end

%This function is used for SFS problem bound chacking 
function V = Bound_Checking(p,lb,ub)
for i=1:size(p,1)
    Flag4ub=p(i,:)>ub;
    Flag4lb=p(i,:)<lb;
    p(i,:)=(p(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
end
V = p;
end

%This function is considered for plotting Positions movement and fitness values
function [] = PlotFunction(New_Position, Lband, Uband,F,G)
    fplot(1:G,1) = sort(F,'descend');
    figure(1)
    clf
    hold on
    %Plotting Positions Movement
    view(10,10)
    set(gca,'XLim',[Lband, Uband],'YLim',[Lband, Uband], 'ZLim',[Lband Uband])
    subplot(2,1,1)
    plot3(New_Position(:,1),New_Position(:,2), New_Position(:,3),'oblack','MarkerFaceColor','g');
    title('Positions Movement')
    grid on;
    %Plotting the value of fitness function
    set(gca,'XLim',[Lband Uband],'YLim',[Lband Uband])
    subplot(2,1,2)
    semilogy(fplot(1:G,1),'--.black')
    title('Fitness Value')
    refreshdata
    drawnow
end
