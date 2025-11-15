function [ BestNeighborPosition ] = CA( iIndex,X,N,dim,fobj,type )
%CA 此处显示有关此函数的摘要
%   此处显示详细说明
BestNeighborFitness=inf;
BestNeighborPosition=zeros(1,dim);
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
            if fobj(X((i-1)*5+ColumnIndex,:))<BestNeighborFitness
                BestNeighborPosition=X((i-1)*5+ColumnIndex,:);
            end
        end
    end
    for j=1:ColumnNum
        if abs(j-ColumnIndex)==1%%邻域左右邻
            if fobj(X((RowIndex-1)*5+j,:))<BestNeighborFitness
                BestNeighborPosition=X((RowIndex-1)*5+j,:);
            end
        end        
    end    
elseif type==2 %%L9
    for i=1:RowNum
        if abs(i-RowIndex)<=2&&abs(i-RowIndex)>0%%邻域上下邻
            if fobj(X((i-1)*5+ColumnIndex,:))<BestNeighborFitness
                BestNeighborPosition=X((i-1)*5+ColumnIndex,:);
            end
        end
    end
    for j=1:ColumnNum
        if abs(j-ColumnIndex)<=2&&abs(j-ColumnIndex)>0%%邻域左右邻
            if fobj(X((RowIndex-1)*5+j,:))<BestNeighborFitness
                BestNeighborPosition=X((RowIndex-1)*5+j,:);
            end
        end        
    end
elseif type==3 %%C9
    for i=1:RowNum
        for j=1:ColumnNum
            if abs(i-RowIndex)<=1 && abs(j-ColumnIndex)<=1
                if i~=RowIndex || j~=ColumnIndex
                    if fobj(X((i-1)*5+j,:))<BestNeighborFitness
                        BestNeighborPosition=X((i-1)*5+j,:);
                    end                    
                end
            end
        end
    end
elseif type==4 %%C13
    for i=1:RowNum
        if abs(i-RowIndex)<=2&&abs(i-RowIndex)>0%%邻域上下邻
            if fobj(X((i-1)*5+ColumnIndex,:))<BestNeighborFitness
                BestNeighborPosition=X((i-1)*5+ColumnIndex,:);
            end
        end
    end
    for j=1:ColumnNum
        if abs(j-ColumnIndex)<=2&&abs(j-ColumnIndex)>0%%邻域左右邻
            if fobj(X((RowIndex-1)*5+j,:))<BestNeighborFitness
                BestNeighborPosition=X((RowIndex-1)*5+j,:);
            end
        end        
    end
    for i=1:RowNum
        for j=1:ColumnNum
            if abs(i-RowIndex)==1 && abs(j-ColumnIndex)==1
                if fobj(X((i-1)*5+j,:))<BestNeighborFitness
                    BestNeighborPosition=X((i-1)*5+j,:);          
                end
            end
        end
    end
elseif type==5 %%C21
    for i=1:RowNum
        for j=1:ColumnNum
            if abs(i-RowIndex)<=2 && abs(j-ColumnIndex)<=2
                if (i~=RowIndex || j~=ColumnIndex) && (abs(i-RowIndex)~=2 || abs(j-ColumnIndex)~=2)
                    if fobj(X((i-1)*5+j,:))<BestNeighborFitness
                        BestNeighborPosition=X((i-1)*5+j,:);
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
                    if fobj(X((i-1)*5+j,:))<BestNeighborFitness
                        BestNeighborPosition=X((i-1)*5+j,:);
                    end                    
                end
            end            
        end
    end
end
end

