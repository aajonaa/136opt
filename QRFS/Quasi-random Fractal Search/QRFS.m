%% General Information
%   关注微信公众号：优化算法侠 Swarm-Opti，获得更多免费、精品、独创、定制代码                                                       %   关注微信公众号：优化算法侠 Swarm-Opti，获得更多免费、精品、独创、定制代码     

% Quasi-random Fractal Search (QRFS)
% This is the code for the QRFS propused in:
% doi: https://doi.org/10.1016/j.eswa.2024.124400
% Authors: Luis A. Beltran, Mario A. Navarro, Diego Oliva, Diego Campos-Peña, Jorge Armando Ramos-Frutos, Saúl Zapotecas-Martínez
% Email: diego.oliva@cucei.udg.mx
% Afer use the code please cite the main paper of QRFS:
% Quasi-random Fractal Search (QRFS): A dynamic metaheuristic with sigmoid population decrement for global optimization
% Expert Systems with Applications
% Volume 254, 15 November 2024, 124400



close all
clear all
%% Set parameters

% Initialize empty instances
Popilations = [];
finalVector = [];

% population size per fractal
npop = 100;                 % Population base size
popIni = npop;              % Initial population size
popFi = 20;                 % Final population size

% Global parameters
dims = 2;                   % Dimensions
nfractals = 5;              % Number of fractals

% limits/boundaries
% These limits were defined for Ackley function
bounds = zeros(dims,2);
bounds(:,1) = -30;
bounds(:,2) = 30;

lb = bounds(:,1);           % lower bound
ub = bounds(:,2);           % upper bound
FName = @ackley;            % Objetive function, here you can change according to the problem to solve

% Number of iterations
maxIter = 1000;

% Variables for initializtion
secuencia = 2;       % Chose a random LDS
max_pop = npop*nfractals;   % General population size considering 

%%  Graph option 
% Fractals Behavior (only two dims)
GraphicParticles=false;

% Final curves
FinalPlots = true;

%% Empty arrays
BestIndIter=zeros(maxIter,dims);        % best solution per iteration
BestFitIter=ones(maxIter,1)*Inf;        % best fit per iteration
FAccessIter=zeros(maxIter,1);           % function accesses per iteration
FAccess=0;                              % function accesses counter
PopulationDecrease = zeros(1,maxIter);  % Population size per iteration

%% Initialization
pop = init_pop(max_pop,lb,ub,dims,secuencia);
SubPoblaciones = mat2cell(pop, repmat(npop, 1, nfractals), dims);

for i=1:nfractals
    % Auxiliar fitness array for each fractal
    Fitness = zeros(npop,1);

    % Evaluation
    for j=1:npop
        PobXfractal = SubPoblaciones{i};
        Popilations(i).Pop = PobXfractal;
        Fitness(j)=feval(FName,PobXfractal(j,:));
        FAccess = FAccess + 1;
    end

    % Best particle (centroid) Selection
    Popilations(i).Fitness = Fitness;
    [MinP, Index_Mim] = min(Popilations(i).Fitness);
    BestFitness(i)  = MinP;
    Popilations(i).MejorMin = MinP;
    [MaxP, Index_Max] = max(Popilations(i).Fitness);
    cen = Popilations(i).Pop(Index_Mim,:);
    Lub = Popilations(i).Pop(Index_Max,:);
    Dis = (Lub - cen);
    Lim1 = cen - Dis;
    Lim2 = cen + Dis;
    % New Bounds
    Lim1 = min(Lim1,ub);
    Lim1 = max(Lim1,lb);
    Lim2 = min(Lim2,ub);
    Lim2 = max(Lim2,lb); 
    
    % update
    Popilations(i).Centroide = cen;
    Popilations(i).LowerBound = Lim2;
    Popilations(i).UpperBound = Lim1;
    finalVector = vertcat(finalVector, Popilations(i).Pop);
end

[BestFit, index] = min(BestFitness);
BestCentoide = Popilations(index).Centroide;

display(['The best position is: ', num2str(BestCentoide)]);
display(['The best fitness is: ', num2str(BestFit)]);

%% Iterative process
for i=1:maxIter
% Calculate population size based on the initial population, iteration
% number and final size population
npop = floor(1+popIni-(1/((1/(popIni-popFi))+exp(-(popFi/maxIter)*i))));

    % 
    for j=1:nfractals
        % new fractal
        Newpop = init_pop(npop,Popilations(j).LowerBound,Popilations(j).UpperBound,dims,randi(4));
        Popilations(j).Pop = Newpop;
        FitnessFractal = zeros(npop,1);

        % Evaluation
        for k=1:npop
            FitnessFractal(k)=feval(FName,Popilations(j).Pop(k,:));
            FAccess = FAccess + 1;
        end

        % Local Selection
        Popilations(j).Fitness = FitnessFractal;
        FitnessFractal = [];
        [MinP, Index_Mim] = min(Popilations(j).Fitness);
        BestFitness(j)  = MinP;
        Popilations(j).MejorMin = MinP;
        cen = Popilations(j).Pop(Index_Mim,:);
        Popilations(j).Centroide = cen;

        % Best particle (centroid) Selection
        NewCentroide = Popilations(j).Centroide + 0.3*(BestCentoide - Popilations(j).Centroide);
        [MaxP, Index_Max] = max(Popilations(j).Fitness);
        Lub = Popilations(j).Pop(Index_Max,:);
        Dis = (Lub - NewCentroide);
        Lim1 = NewCentroide - Dis;
        Lim2 = NewCentroide + Dis;
        % New Bounds
        Lim1 = min(Lim1,ub);
        Lim1 = max(Lim1,lb);
        Lim2 = min(Lim2,ub);
        Lim2 = max(Lim2,lb); 

        % Update
        Popilations(j).Centroide = NewCentroide;
        Popilations(j).LowerBound = Lim2;
        Popilations(j).UpperBound = Lim1;
        finalVector = vertcat(finalVector, Popilations(j).Pop);
    end

    % Global Selection
    [GlobalFractalFitness(i), index] = min(BestFitness);
    BestCentoide = Popilations(index).Centroide; % Best solution

    % for j=1:nfractals
    %     NewCentroide = Poblaciones(j).Centroide + 0.3*(BestCentoide - Poblaciones(j).Centroide);
    %     [MaxP, Index_Max] = max(Poblaciones(j).Fitness);
    %     Lub = Poblaciones(j).Pop(Index_Max,:);
    %     Distancia = (Lub - NewCentroide);%/2;
    %     Limite1 = NewCentroide - Distancia;
    %     Limite2 = NewCentroide + Distancia;
    %     %Limites
    %     Limite1 = min(Limite1,ub);
    %     Limite1 = max(Limite1,lb);
    %     Limite2 = min(Limite2,ub);
    %     Limite2 = max(Limite2,lb); 
    %     %%%%%%%%%%%%%%%%%%%
    %     Poblaciones(j).Centroide = NewCentroide;
    %     Poblaciones(j).LowerBound = Limite2;
    %     Poblaciones(j).UpperBound = Limite1;
    %     vectorFinal = vertcat(vectorFinal, Poblaciones(j).Pop);
    % end
   PopulationDecrease(i) = npop;
   
   [bestever,~]=min(BestFitness);
   
    %% GraphicParticles
    if (GraphicParticles==true)
        %Plot
        for h=1:npop
            resolution=((ub(1)-lb(1))/2)/20;
            [X,Y] = meshgrid(lb(1):resolution:ub(1),lb(2):resolution:ub(2));
            %colormap(winter);
            CM=[0.3:1/20:1]';
            CM=[CM,CM,CM];
            colormap(CM);
            for k=1:size(X,1)
                for l=1:size(X,2)
                    Z(k,l) = feval(FName,[X(k,l),Y(k,l)]);
                end
            end
            a= gradient(Z);
            surf(X,Y,Z,'AlphaData',a,'FaceAlpha',.3,'EdgeAlpha',.5);
            axis([lb(1) ub(1) lb(2) ub(2)]);
            hold on
            for nf=1:nfractals
                plot3(Popilations(nf).Pop(1),Popilations(nf).Pop(2),Fitness(h),'Marker','o','MarkerFaceColor',[0 .8 1],'MarkerSize',4,'LineWidth',.5,'MarkerEdgeColor','k');
            end
            hold off;
            view(0,-90)
            pause(0.01);
        end
    end
   
    % update
    BestFitIter(i)= bestever;
    BestIndIter(i,:)=BestCentoide;
    FAccessIter(i) = FAccess;
end

% Best solution variables
bestFit = bestever;
bestInd = BestCentoide;


% Plots
if FinalPlots
    figure()
    semilogy(BestFitIter)
    title('Convergence curve semylogy') % Convergence curve
    
    figure()
    loglog(GlobalFractalFitness) 
    title('Best fitness per iteration loglog') % Best fitness per iteration
    
    figure()
    plot(PopulationDecrease) 
    title('Population decreasement') % Population decreasement
end