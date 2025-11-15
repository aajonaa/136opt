%**************************************************************************
% Stochastic Fractal Search Algorithm                                     %
% The output function is:                                                 %
% pbest: the best solution for problem                                    %
% fbest: the value of best fitness function                               %
% F: recording the value of best fitness function in each iteration       %
%**************************************************************************
function [best_pos,best_score,Convergence_curve] = SFS(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
%Creating random Positionss in considered search space=========================
% Positions = repmat(S.Lband,S.Start_Positions,1) + rand(S.Start_Positions, S.Ndim).* ...
%    (repmat(S.Uband - S.Lband,S.Start_Positions,1));
best_pos = zeros(1,dim);
best_score = inf; %change this to -inf for maximization problems
%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
%==========================================================================
%Calculating the fitness of first created Positionss===========================
FitnessVal = [];
for i = 1 : size(Positions,1)
    FitnessVal(i,:) = fobj(Positions(i,:));
end
[Sorted_FitnessVal, index] = sort(FitnessVal);
Positions = Positions(index,:);%sorting the Positionss based on obtaind result
%==========================================================================
%Finding the Best Positions in the group=======================================
best_pos = Positions(1, :);
best_score = Sorted_FitnessVal(1);%saving the first best fitness
F = Sorted_FitnessVal(1);
%==========================================================================
%Starting Optimizer========================================================
FEs = 0;
it = 1;
Maximum_Diffusion = 2;
Walk = 0.75;%也可以为 0.75
%if you want to see the best result in each iteration, set it to 1.
ShowResult = 0;
while FEs < MaxFEs
    %for G = 1  : S.Maximum_Generation
    New_Positions = [];%creating new Positions
    FitVector = [];%creating vector of fitness functions
    %diffusion process occurs for all Positionss in the group
    for i = 1 : size(Positions,1)
        %creating new Positionss based on diffusion process
        [NP, fit] = Diffusion_Process(Positions(i,:),Maximum_Diffusion,Walk,lb,ub,FEs,best_pos,fobj);
        New_Positions = [New_Positions;NP];
        FitVector = [FitVector,fit];
    end
    %======================================================================
    
    %updating best Positions obtained by diffusion process
    best_score = min(FitVector);
%     if temp_best_score < best_score
%         best_score = temp_best_score;
%         best_pos = New_Positions(find(FitVector == best_score),:);
%     end
    best_pos = New_Positions(find(FitVector == best_score),:);
    SearchAgents_no = size(New_Positions,1);
    fit = FitVector';
    [sortVal, sortIndex] = sort(fit);
    
    %Starting The First Updating Process====================================
    for i=1:SearchAgents_no
        Pa(sortIndex(i)) = (SearchAgents_no - i + 1) / SearchAgents_no;
    end
    RandVec1 = randperm(SearchAgents_no);
    RandVec2 = randperm(SearchAgents_no);
    
    for i = 1 : SearchAgents_no
        for j = 1 : size(New_Positions,2)
            if rand > Pa(i)
                P(i,j) = New_Positions(RandVec1(i),j) - rand*(New_Positions(RandVec2(i),j) - ...
                    New_Positions(i,j));
            else
                P(i,j)= New_Positions(i,j);
            end
        end
    end
    P = Bound_Checking(P,lb,ub);%for checking bounds
    Fit_FirstProcess = [];
    for i = 1 : SearchAgents_no
        FEs = FEs + 1;
        Fit_FirstProcess = [Fit_FirstProcess;fobj(P(i,:))];
    end
    for i=1:SearchAgents_no
        if Fit_FirstProcess(i,:) <= fit(i,:)
            New_Positions(i,:)=P(i,:);
            fit(i,:)=Fit_FirstProcess(i,:);
        end
    end
    FitVector = fit;
    %======================================================================
    
    [SortedFit,SortedIndex] = sort(FitVector);
    New_Positions = New_Positions(SortedIndex,:);
    best_pos = New_Positions(1,:);%first Positions is the best
    temp1_best_score = SortedFit(1,1); %%
    if temp1_best_score < best_score
        best_score = temp1_best_score;
        best_pos = New_Positions(1,:);%first Positions is the best
    end
    F = [F;FitVector(1,1)];
    F = sort(F);
    %Convergence_curve(it)=best_score;
    %it = it + 1;
    
    %temp1_best_score = FitVector(SortedIndex(1),1);
    %if temp1_best_score < best_score
        %best_score = temp1_best_score;
    %end
%     Convergence_curve(it) = best_score;
%     it = it + 1;
%    if ShowResult == 1
%         fprintf('Iteration: %i', G);
%         fprintf(',    Best result: %e \n', F(1,1));
%     end
    
    %Plotting All Positionss===================================================
    %if eq(S.plot,1)
    %    PlotFunction(New_Positions, min(S.Lband), max(S.Uband),F,G+1);
    %end
    %======================================================================
    
%     best_score = SortedFit(1,:);
%     if best_score <= F(1,:)
%         best_pos = New_Positions(1,:);
%         best_score = F(1,:);
%     end
    
    Positions = New_Positions;
    
    %Starting The Second Updating Process==================================
    Pa = sort(SortedIndex/SearchAgents_no, 'descend');
    
    for i = 1 : SearchAgents_no
        if rand > Pa(i)
            %selecting two different Positionss in the group
            R1 = ceil(rand*size(Positions,1));
            R2 = ceil(rand*size(Positions,1));
            while R1 == R2
                R2 = ceil(rand*size(Positions,1));
            end
            
            if rand < .5
                ReplacePositions = Positions(i,:) - ...
                    rand * (Positions(R2,:) - best_pos);
                ReplacePositions = Bound_Checking(ReplacePositions,lb,ub);
            else
                ReplacePositions = Positions(i,:) + ...
                    rand * (Positions(R2,:) - Positions(R1,:));
                ReplacePositions = Bound_Checking(ReplacePositions,lb,ub);
            end
            FEs = FEs + 1;
            if fobj(ReplacePositions) < fobj(Positions(i,:))
                Positions(i,:) = ReplacePositions;
            end
            %if feval(S.Function_Name, ReplacePositions) < ...
            %        feval(S.Function_Name, Positions(i,:))
            %    Positions(i,:) = ReplacePositions;
        end
    end
    %best_score = min(F);
    Convergence_curve(it) = best_score;
   % display(best_score);
    it = it + 1;
end
%======================================================================
end

%This function is used for SFS problem bound chacking
function V = Bound_Checking(p,lb,ub)
%     for i = 1 : size(p,1)
%         upper = double(gt(p(i,:),upB));
%         lower = double(lt(p(i,:),lowB));
%         up = find(upper == 1);
%         lo = find(lower == 1);
%         if (size(up,2)+ size(lo,2) > 0 )
%             for j = 1 : size(up,2)
%                 p(i, up(j)) = (upB(up(j)) - lowB(up(j)))*rand()...
%                     + lowB(up(j));
%             end
%             for j = 1 : size(lo,2)
%                 p(i, lo(j)) = (upB(lo(j)) - lowB(lo(j)))*rand()...
%                     + lowB(lo(j));
%             end
%         end
%     end
for i=1:size(p,1)
    Flag4ub=p(i,:)>ub;
    Flag4lb=p(i,:)<lb;
    p(i,:)=(p(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
end
V = p;
end