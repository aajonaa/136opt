function colormap = createnewColormap(colors, num)
%  --------------------------------------------------------
% CREATECOLORMAP Create a new colormap based on the specified,
%                                 color vectors and numbers
%
%   INPUT:
%      colors:  color vectors
%      num:    color numbers
%
%   OUTPUT:     
%      colormap: a new colormap
% --------------------------------------------------------
m = size(colors, 1);
if m >= num
    colormap = colors;
else
    % Reset range
    range = 0 : m-1;
    range = range*(num-1)/(m-1) + 1;
    % Interpolate
    colormap = nan(num, 3);
    for i = 1:3
        colormap(:, i) = interp1(range, colors(:, i), 1:num);
    end
end