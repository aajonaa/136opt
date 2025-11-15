cons = table(connections);
writetable(cons, 'connections3.xlsx');

% 创建连接矩阵转换脚本

% 以下是示例矩阵，您可以替换为您的实际矩阵
% connectionMatrix = rand(30,30);
% connectionMatrix = round(connectionMatrix);  % 转换为0-1矩阵
connectionMatrix = readtable("connections3.xlsx");
connectionMatrix = table2array(connectionMatrix);

% 确保矩阵是对称的并且对角线为0
connectionMatrix = triu(connectionMatrix, 1);

% 初始化一个空的cell数组来存储连接
nodeConnections = {};

% 遍历上三角矩阵（不包括对角线）
for i = 1:size(connectionMatrix, 1)-1
    for j = i+1:size(connectionMatrix, 2)
        if connectionMatrix(i, j) == 1
            % 添加连接关系到cell数组
            nodeConnections{end+1, 1} = ['Node' num2str(i)];
            nodeConnections{end, 2} = ['Node' num2str(j)];
        end
    end
end

% 将connections转换为表格
connectionsTable = cell2table(nodeConnections, 'VariableNames', {'Source', 'Target'});

% 保存到Excel
writetable(connectionsTable, 'ConnectionMatrixBased3.xlsx');

% 在命令窗口显示连接关系
disp('连接关系:');
disp(connectionsTable);

% 打印连接数
fprintf('总连接数: %d\n', size(nodeConnections, 1));