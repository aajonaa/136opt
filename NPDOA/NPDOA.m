%------------------------------- Reference --------------------------------
% T. X. Wu, A brain-inspired neuron population dynamic optimization algorithm. MIT Press,
% 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

function [best_pos,convergence_curve] = NPDOA(N,MaxFEs,Lb,Ub,dim,fobj)

    %% Parameter setting
    AT = 0.1;
    Nd = 0.1;
    a = 0.7;
    d = 0.1;

    FEs = 0;
    
    %% Generate random population
    X = initialization(N,dim,ub,lb);
    it = 1;

    AllFitness = inf * ones(1, N);

    for i = 1:N
        AllFitness(1, i) = fobj(X(i, :));
        FEs = FEs + 1;
    end

    [~, ind] = min(AllFitness);
    bestFitness = min(AllFitness);
    best_pos = X(ind, :);
    
    %% Optimization
    while FEs < MaxFEs
        % 将最优个体视为吸引子，将附近的神经状态活动模式拉向稳定状态

        % 将多个最优个体视为吸引子集，将附近的神经状态活动模式不同程度的吸引
        [~,rank]   = sort(AllFitness(X),"ascend");
        X = X(rank);
%                 [~,best]   = min(FitnessSingle(Population));
%                 Gbest      = Population(best);
        attractordec = X(1:floor(end*AT));
        % 由神经群网络组成的动力学模型通常假设神经群之间的相互作用是通过加性耦合或扩散耦合进行的。
        % 当使用加性耦合时，种群的活动受到邻近种群活动总和的影响。相反，当使用扩散耦合时，神经种群受到其活动与其邻居活动之差之和的影响
        Offspring = X.decs;
        T = Problem.FE/Problem.maxFE;
        for i = floor(Problem.N*AT):Problem.N
            newone = X(i).decs;
            for j = 1:floor(Problem.N*AT)
            % 随机创建该神经群与吸引子的邻接矩阵
%                         A = rand(1,Problem.D);
%                         A((A<T))=1;
%                         A((A~=1))=0;
                A = randi([0,1] ,1, Problem.D);
                SA = randi([1, floor(Problem.N*AT)]);
                bestdec = attractordec(SA).decs;
%                     A=randi([1,1],1,Problem.D);
                newdec = X(i).decs;
                L = rand(1, Problem.D);
                L1= rand(1, Problem.D);
                newone = newone + a .*L.* L .*A .* (bestdec-newdec) ;
%                     newone = newdec;
            end
            % 高斯噪声
            Lower = min(X.decs);
            Upper = max(X.decs);
            noiseMean = 0 * ones(1,Problem.D);
            noiseStd = Nd*(Upper - Lower);
%                     noiseStd = Nd;
%                     mu = rand(1,1);
%                     if mu > 0.99
%                         noiseStd = noiseStd * 2;
%                     end
%                     noiseStd = noiseStd * (1-Problem.FE/Problem.maxFE) + 0.1;
            noise = randn() * noiseStd + noiseMean;
            newone = newone + noise;
            % 随机创建神经群间的接矩阵
%                     Newdecs = ones(Problem.N, Problem.D);
%                     for j = 1:Problem.N
%                         Newdecs(j,:) = newdec;
%                     end
             Offspring(i,:) = newone;
        end

        meanN = mean(Offspring);
        for i = 1:Problem.N
%                     addcontact = zeros(1, Problem.D);
%                     difcontact = zeros(1, Problem.D);
            newone = Offspring(i,:);
%                     for j = 1:Problem.N
%                         tempdec = Offspring(j,:);
%                         A = randi([0,1] ,1, Problem.D);
% %                         A = randi([1,1],1,Problem.D);
%                         tempadd = L.* A.*(tempdec-newone);
% %                         tempdif =  A.*(tempdec);
% % 
%                         addcontact = addcontact + tempadd;
% %                         difcontact = difcontact + tempdif;
%                     end
            
            l1 = 1/2*d*rand(1, Problem.D);
            l2 = 1/2*d*rand(1, Problem.D);
            L= randi([-1,1] ,1, Problem.D);
            L1= randi([-1,1] ,1, Problem.D);
            L2= randi([-1,1] ,1, Problem.D);
            A = randi([0,1] ,1, Problem.D);
            A1 = randi([0,1] ,1, Problem.D);
            A2 = randi([0,1] ,1, Problem.D);
%                     A1=A2;
%                     L1 = L2;
%                     A = rand(1,Problem.D);
%                     A((A<T))=1;
%                     A((A~=0))=1;
%                     newone = newone + L.*addcontact/Problem.N + L.*difcontact/Problem.N;
%                     newone = newone + L.*addcontact/Problem.N;
            r = rand();
%                     if r >= 0.5
            
            newone = newone + (L1.*l1.* A1.* (newone-meanN)+L2.*A2.*l2.*meanN) * rand()*(1-T);

%                     newone = newone + (L1.*l1.* A.* (newone-meanN)) * rand()*(1-T);
%                     else
%                         newone = newone - L1.* L.*A.* (newone-meanN)* (1-T)*rand();
%                     end
            
            Offspring(i,:) = newone;
        end

        Offspring  = Problem.Evaluation(Offspring);

        X = Offspring;
    end
end
