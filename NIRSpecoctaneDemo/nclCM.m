function [colorList,colorName]=nclCM(varargin)
% @author: slandarer/Zhaoxu Liu
% The NCAR Command Language LIB (NCL), NCL is a scripting language
% for the analysis and visualization of climate and weather data.
%
% The default parameter is to return a list of all 487 colormap identifiers
% Data source from:
% Zhaoxu Liu / slandarer (2022). ncl colormaps 
% https://www.mathworks.com/matlabcentral/fileexchange/120848-ncl-colormaps, 
% MATLAB Central File Exchange,2022/11/19.
% -----------------------
%   INPUT:
%      type:  Color type can be string or index number,490 in all
%      num:  Color scale or level, default is 256
%
%   OUTPUT:
%      Color: A specified color matrix (double)
%      colorName: A specified color name (string)
%
%  Example:
%	   nclCM; %show colorbar image of all color
%      [Color,colorName] = nclCM(60);
% 	   [Color,colorName] = nclCM(60,256);
% 	   [Color,colorName] = nclCM(60,20); 
% --------------------------------------------------------
narginchk(0,2)

nclCM_Data = load('.\NIRSpecoctaneDemo\colorData\nclCM_Data.mat');
CList_Data   = nclCM_Data.Colors;

if nargin<2
     num=-1;
else
     num = varargin{2};
end

% If the input parameter is empty, execute the display of all color bars
if nargin<1     
    colorName = nclCM_Data.Names;
    %disp(ctype);
    dispallBars(colorName,411,490,'ncl'); % show all bars
    return
    % error('nclcm:missingarguments',...
    % 'Specify the type or index of colormap');
end

type = varargin{1};
if isnumeric(type)
    if type>490
        error('dipcm:agarguments','The color index range is [1~490]');
    end
    Cmap = CList_Data{type};
    colorName =nclCM_Data.Names{type};
else
    Cpos  = strcmpi(type,nclCM_Data.Names);
    Cmap = CList_Data{find(Cpos,1)};
    colorName = type;
end

% Interpolate to the given number of colors
if num>0
    colorList = interpCMap(Cmap, num);
else
    colorList = Cmap;
end
end