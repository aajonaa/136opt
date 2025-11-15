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

function pop = init_pop(pop_size, varmin, varmax, dim, secc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input values:
%pop_size = number of solutions/particles/individuals
%varmin = lower limit
%varmax = upper limit
%dim = number of dimensions
%secc = Low discrepancy sequence selected

% To select the desired sequence 
% it is necessary to choose a number based on:
% 1.- Sobol
% 2.- Latin Hypercube
% 3.- Halton
% 4.- Hammerley
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    VarMax = varmax;
    VarMin = varmin;
    nPop = pop_size;
    dims = dim;
    
    % pop = randn(10,2).*10;
    empty_individual.Position=[];
    empty_individual.Cost=[];
    
    pop=repmat(empty_individual,nPop,1);
    
    for i=1:nPop
        for ard=1:dims
            pop(i).Position(ard)=unifrnd(VarMin(ard),VarMax(ard));
        end
    end
    
    %sobol
    if secc == 4
        pop_s = SobolInicialization(dims, VarMin, VarMax, nPop,pop);
        pop = vertcat( pop_s.Position );
    end
    
    %latin
    if secc == 2
        pop_l = LatinHyperCubeInicialization(dims, VarMin, VarMax, nPop,pop);
        pop = vertcat( pop_l.Position );
    end
    
    %Halton
    if secc == 3
        pop_h = HaltonInicialization(dims, VarMin, VarMax, nPop,pop);
        pop = vertcat( pop_h.Position );
    end
    
    % %Hammer
    if secc == 1
        pop_hm = HammerleyInicialization(dims, VarMin, VarMax, nPop,pop);
        pop = vertcat( pop_hm.Position );
    end

end