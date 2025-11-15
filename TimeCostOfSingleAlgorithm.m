clear all 
close all
clc

algorithmName = {'ALCPSO', 'SBO'};

% 'CEC2017'
Function_name_all={'F107','F108','F109','F110','F111','F112','F113','F114','F115','F116','F117','F118','F119','F120','F121','F122','F123','F124','F125','F126','F127','F128','F129','F130','F131','F132','F133','F134','F135','F136'};

for dim = [10 30]
    dim
    for algo = 1:size(algorithmName, 2)
        algo_fhd = str2func(algorithmName{algo})
    
        tic
        for funcNum=1:size(Function_name_all,2)
            if funcNum == 2
                continue;
            end
        
            fhd = @(x) cec17_func(x', funcNum);
            parfor cflod=1:4
                % display(['flod',num2str(cflod)]);
                [best_pos,cg_curve]=algo_fhd(30,300000,-100,100,dim,fhd);
            end
        end
        toc
        Time = toc;
    end
end
