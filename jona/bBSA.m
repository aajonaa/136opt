%{
Backtracking Search Optimization Algorithm (BSA)

Platform: Matlab 2013a   


Cite this algorithm as;
[1]  P. Civicioglu, "Backtracking Search Optimization Algorithm for 
numerical optimization problems", Applied Mathematics and Computation, 219, 8121?144, 2013.


Copyright Notice
Copyright (c) 2012, Pinar Civicioglu
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.


%}

function  [best_score, best_pos,CNVG, Time] = bBSA(N, epoch, dim, A, trn, vald, ~, classifierFhd)

tic;
% Maybe there is no need to copy this code.
if (nargin<8)
    str = 'knn';
    classifierFhd = Get_Classifiers(str);
end

TFid = 1;

DIM_RATE = 1;  %此变量的取值在原文中并没有给出
CNVG=[];
FES = 0;
% Should the best score be inf.
% best_score = 1;
best_score = inf;
best_pos = zeros(1,dim);
it = 1;
fitnessX = inf * ones(N, 1);
%INITIALIZATION
% if numel(low)==1, low=low*ones(1,dim); up=up*ones(1,dim); end % this line must be adapted to your problem
%numel:Number of elements in an array or subscripted array expression.
X=initialization(N,dim,1,0) > 0.5; % see Eq.1 in [1]
% fitnessX=feval(fobj,X);
for i = 1 : N
    % fitnessX(i) = feval(fobj,X(i,:));
    fitnessX(i) = AccSz2(X(i,:), A,trn,vald,classifierFhd);
    FES  = FES +1;
    if fitnessX(i) < best_score
        best_score = fitnessX(i);
        best_pos = X(i,:);
    end
end
historical_X = initialization(N, dim, 1, 0) > 0.5;% see Eq.2 in [1]

% historical_X  is swarm-memory of BSA as mentioned in [1].

% ------------------------------------------------------------------------------------------ 
while FES < epoch
    %SELECTION-I
    if rand<rand, historical_X=X; end  % see Eq.3 in [1] redefine oldP at the beginning of each iteration
    historical_X=historical_X(randperm(N),:); % see Eq.4 in [1] randomly change the order of the individuals in oldP
    F=get_scale_factor; % see Eq.5 in [1], you can define other F generation strategies 
    map=zeros(N,dim); % see Algorithm-2 in [1]  我看论文中这里是定义的全为1的矩阵呢？       
    if rand<rand,
        for i=1:N,  u=randperm(dim); map(i,u(1:ceil(DIM_RATE*rand*dim)))=1; end
    else
        for i=1:N,  map(i,randi(dim))=1; end
    end
    % RECOMBINATION (MUTATION+CROSSOVER)   
    offsprings=X+(map.*F).*(historical_X-X);   % see Eq.5 in [1]  
    %原文中不是X+F*(historical_X-X)吗？对，但这里将变异与交叉合并在一起了，  
    % offsprings=BoundaryControl(offsprings,low,up); % see Algorithm-3 in [1]
    for i = 1:N
        for j = 1:dim
            offsprings(i, j) = trnasferFun(offspring(i, j), offsprings(i, j), TFid);
        end
    end

    % SELECTON-II
    for i = 1:N
        % fitnessoffsprings(i) = feval(fobj,offsprings(i,:));
        fitnessoffsprings(i) = AccSz2(offsprings(i,:), A,trn,vald,classifierFhd);
        FES = FES + 1;
        if fitnessoffsprings(i) < best_score
            best_score = fitnessoffsprings(i);
            best_pos = offsprings(i,:);
        end
    end
    ind=fitnessoffsprings<fitnessX;
    fitnessX(ind)=fitnessoffsprings(ind);%讲fitnessoffspring中适应度值小的赋值给fitnessX中相应位置
    X(ind,:)=offsprings(ind,:);
    [globalminimum,ind]=min(fitnessX);    
    globalminimizer=X(ind,:);
    best_score = globalminimum;
    best_pos =globalminimizer;
    CNVG(it)=best_score;
    it =it + 1;
end

Time = toc;

end

% function X=BoundaryControl(X,low,up)
% [N,dim]=size(X);
% for i=1:N
%     for j=1:dim                
%         k=rand<rand; % you can change boundary-control strategy
%         if X(i,j)<low(j)
%             if k, X(i,j)=low(j); 
%             else X(i,j)=rand*(up(j)-low(j))+low(j); 
%             end 
%         end        
%         if X(i,j)>up(j)
%             if k, X(i,j)=up(j);  
%             else
%                 X(i,j)=rand*(up(j)-low(j))+low(j); 
%             end 
%         end
%     end
% end
% % return
% end



function F=get_scale_factor % you can change generation strategy of scale-factor,F    
     F=3*randn; % STANDARD brownian-walk
    % F=4*randg;  % brownian-walk    
    % F=lognrnd(rand,5*rand);  % brownian-walk              
    % F=1/normrnd(0,5);        % pseudo-stable walk (levy-like)
    % F=1./gamrnd(1,0.5);      % pseudo-stable walk (levy-like, simulates inverse gamma distribution; levy-distiribution)   
% return

end