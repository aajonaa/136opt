% This function opens multiple .fig files sequentially
% based on the number of functions specified in funcNum.

funcNum = 29;

for func = 1:funcNum
    if func == 2
        continue
    end
    % Create the filename for the .fig file
    figFile = ['SBO-F', num2str(func), '-CC.fig'];
    
    % Open the .fig file
    openfig(figFile);
    
    % Optional: display a message to confirm the file is opened
    fprintf('Opened file: %s\n', figFile);
end
