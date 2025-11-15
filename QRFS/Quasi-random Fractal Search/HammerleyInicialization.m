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

 function [pop]  = HammerleyInicialization(dim,down,up,nPop,pop)
nPopBase=nPop;
if (nPopBase>100)
    nPopBase=100;
else
    nPopBase=nPop;
end 
H = hammersley_sequence ( 1, nPop, dim, nPopBase );
H=H';
% scatter(X0(:,1),X0(:,2),20,'b')
% axis square
% title('{\bf Quasi-Random Scatter}')
for i=1:nPop

     for j=1:dim
       pop(i).Position([j])=H(i,j).*(up(j)-down(j))+down(j);
     end

end    
   
 end