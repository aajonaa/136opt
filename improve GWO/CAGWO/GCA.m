function [ Alpha_pos,Beta_pos,Delta_pos ,Fes] = GCA( iIndex,X,N,dim,fobj,type,Fes )%i：行数（第几个灰狼）Positions：灰狼位置信息 SearchAgents_no：灰狼个数 dim：维度 fobj：目标函数  type：类型
%CA 此处显示有关此函数的摘要
%   此处显示详细说明
% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

RowNum=N/5;
ColumnNum=5;
RowIndex=ceil(iIndex/5);
ColumnIndex=rem(iIndex,5);
if ColumnIndex==0
    ColumnIndex=5;
end
if type==1 %%L5
    for i=1:RowNum
        if abs(i-RowIndex)==1%%邻域上下邻
            % Update Alpha, Beta, and Delta
			Fes = Fes + 1;
            NeighborFitness=fobj(X((i-1)*5+ColumnIndex,:));
            if NeighborFitness<Alpha_score 
                Alpha_score=NeighborFitness; % Update alpha
                Alpha_pos=X((i-1)*5+ColumnIndex,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                Beta_score=NeighborFitness; % Update beta
                Beta_pos=X((i-1)*5+ColumnIndex,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                Delta_score=NeighborFitness; % Update delta
                Delta_pos=X((i-1)*5+ColumnIndex,:);
            end
        end
    end
    for j=1:ColumnNum
        if abs(j-ColumnIndex)==1%%邻域左右邻
			Fes = Fes + 1;
            NeighborFitness=fobj(X((RowIndex-1)*5+j,:));
                    % Update Alpha, Beta, and Delta
            if NeighborFitness<Alpha_score 
                Alpha_score=NeighborFitness; % Update alpha
                Alpha_pos=X((RowIndex-1)*5+j,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                Beta_score=NeighborFitness; % Update beta
                Beta_pos=X((RowIndex-1)*5+j,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                Delta_score=NeighborFitness; % Update delta
                Delta_pos=X((RowIndex-1)*5+j,:);
            end
        end        
    end    
elseif type==2 %%L9
    for i=1:RowNum
        if abs(i-RowIndex)<=2&&abs(i-RowIndex)>0%%邻域上下邻
			Fes = Fes + 1;
            NeighborFitness=fobj(X((i-1)*5+ColumnIndex,:));
                    % Update Alpha, Beta, and Delta
            if NeighborFitness<Alpha_score 
                Alpha_score=NeighborFitness; % Update alpha
                Alpha_pos=X((i-1)*5+ColumnIndex,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                Beta_score=NeighborFitness; % Update beta
                Beta_pos=X((i-1)*5+ColumnIndex,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                Delta_score=NeighborFitness; % Update delta
                Delta_pos=X((i-1)*5+ColumnIndex,:);
            end
        end
    end
    for j=1:ColumnNum
        if abs(j-ColumnIndex)<=2&&abs(j-ColumnIndex)>0%%邻域左右邻
			Fes = Fes + 1;
            NeighborFitness=fobj(X((RowIndex-1)*5+j,:));
            % Update Alpha, Beta, and Delta
            if NeighborFitness<Alpha_score 
                Alpha_score=NeighborFitness; % Update alpha
                Alpha_pos=X((RowIndex-1)*5+j,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                Beta_score=NeighborFitness; % Update beta
                Beta_pos=X((RowIndex-1)*5+j,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                Delta_score=NeighborFitness; % Update delta
                Delta_pos=X((RowIndex-1)*5+j,:);
            end
        end        
    end
elseif type==3 %%C9
    for i=1:RowNum
        for j=1:ColumnNum
            if abs(i-RowIndex)<=1 && abs(j-ColumnIndex)<=1
                if i~=RowIndex || j~=ColumnIndex
					Fes = Fes + 1;
                    NeighborFitness=fobj(X((i-1)*5+j,:));
                            % Update Alpha, Beta, and Delta
                    if NeighborFitness<Alpha_score 
                        Alpha_score=NeighborFitness; % Update alpha
                        Alpha_pos=X((i-1)*5+j,:);
                    end

                    if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                        Beta_score=NeighborFitness; % Update beta
                        Beta_pos=X((i-1)*5+j,:);
                    end

                    if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                        Delta_score=NeighborFitness; % Update delta
                        Delta_pos=X((i-1)*5+j,:);
                    end
                end
            end
        end
    end
elseif type==4 %%C13
    for i=1:RowNum
        if abs(i-RowIndex)<=2&&abs(i-RowIndex)>0%%邻域上下邻
			Fes = Fes + 1;
            NeighborFitness=fobj(X((i-1)*5+ColumnIndex,:));
                    % Update Alpha, Beta, and Delta
            if NeighborFitness<Alpha_score 
                Alpha_score=NeighborFitness; % Update alpha
                Alpha_pos=X((i-1)*5+ColumnIndex,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                Beta_score=NeighborFitness; % Update beta
                Beta_pos=X((i-1)*5+ColumnIndex,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                Delta_score=NeighborFitness; % Update delta
                Delta_pos=X((i-1)*5+ColumnIndex,:);
            end
        end
    end
    for j=1:ColumnNum
        if abs(j-ColumnIndex)<=2&&abs(j-ColumnIndex)>0%%邻域左右邻
			Fes = Fes + 1;
            NeighborFitness=fobj(X((RowIndex-1)*5+j,:));
            if NeighborFitness<Alpha_score 
                Alpha_score=NeighborFitness; % Update alpha
                Alpha_pos=X((RowIndex-1)*5+j,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                Beta_score=NeighborFitness; % Update beta
                Beta_pos=X((RowIndex-1)*5+j,:);
            end

            if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                Delta_score=NeighborFitness; % Update delta
                Delta_pos=X((RowIndex-1)*5+j,:);
            end
        end        
    end
    for i=1:RowNum
        for j=1:ColumnNum
            if abs(i-RowIndex)==1 && abs(j-ColumnIndex)==1
				Fes = Fes + 1;
                NeighborFitness=fobj(X((i-1)*5+j,:));
                if NeighborFitness<Alpha_score 
                    Alpha_score=NeighborFitness; % Update alpha
                    Alpha_pos=X((i-1)*5+j,:); 
                end

                if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                    Beta_score=NeighborFitness; % Update beta
                    Beta_pos=X((i-1)*5+j,:); 
                end

                if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                    Delta_score=NeighborFitness; % Update delta
                    Delta_pos=X((i-1)*5+j,:); 
                end
            end
        end
    end
elseif type==5 %%C21
    for i=1:RowNum
        for j=1:ColumnNum
            if abs(i-RowIndex)<=2 && abs(j-ColumnIndex)<=2
                if (i~=RowIndex || j~=ColumnIndex) && (abs(i-RowIndex)~=2 || abs(j-ColumnIndex)~=2)
					Fes = Fes + 1;
                    NeighborFitness=fobj(X((i-1)*5+j,:));
                            % Update Alpha, Beta, and Delta
                    if NeighborFitness<Alpha_score 
                        Alpha_score=NeighborFitness; % Update alpha
                        Alpha_pos=X((i-1)*5+j,:);
                    end

                    if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                        Beta_score=NeighborFitness; % Update beta
                        Beta_pos=X((i-1)*5+j,:);
                    end

                    if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                        Delta_score=NeighborFitness; % Update delta
                        Delta_pos=X((i-1)*5+j,:);
                    end
                end
            end
        end
    end    
else %%C25
    for i=1:RowNum
        for j=1:ColumnNum
            if abs(i-RowIndex)<=2 && abs(j-ColumnIndex)<=2
                if (i~=RowIndex || j~=ColumnIndex)
					Fes = Fes + 1;
                    NeighborFitness=fobj(X((i-1)*5+j,:));
                            % Update Alpha, Beta, and Delta
                    if NeighborFitness<Alpha_score 
                        Alpha_score=NeighborFitness; % Update alpha
                        Alpha_pos=X((i-1)*5+j,:);
                    end

                    if NeighborFitness>Alpha_score && NeighborFitness<Beta_score 
                        Beta_score=NeighborFitness; % Update beta
                        Beta_pos=X((i-1)*5+j,:);
                    end

                    if NeighborFitness>Alpha_score && NeighborFitness>Beta_score && NeighborFitness<Delta_score 
                        Delta_score=NeighborFitness; % Update delta
                        Delta_pos=X((i-1)*5+j,:);
                    end
                end
            end            
        end
    end
end
end

