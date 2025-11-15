function  [Color,colorName] = idipCM(varargin)
%  --------------------------------------------------------
% IDIPCM Define Colormaps for Plots
%
%   The input parameter is marked as cm, and the function returns 
%   256*3 RGB numerical matrix of colormap
%   The default parameter is to return a list of all 80 colormap identifiers
% --------------------------------------------------------
%   INPUT:
%      Type:  Colormap type can be string or index number
%      num: Color scale or level, default is 256
%
%   OUTPUT:
%      Color: A specified color matrix (double)
%      colorName: A specified color name (string)
%  Example:
%	   idipCM; % default as showing colorbar image of all colors
%      [Color,colorName] = idipCM(60); % A specified color with 256 level
% 	   [Color,colorName] = idipCM(60,256); % A specified color with 256
% 	   level
% 	   [Color,colorName] = idipCM(60,20); % A specified color with 20 level
% --------------------------------------------------------
narginchk(0,2)

dipCM_Data = load ('.\NIRSpecoctaneDemo\colorData\dipCM_data.mat');
ColorNames  = dipCM_Data.ColorNames;
ColorsList = dipCM_Data.ColorsList;
%disp(dipCM_Data.Authors);

if nargin >0
    Type = varargin{1};
    if nargin==1
        num = 256;
    else
        num = varargin{2};
    end
    
    if isnumeric(Type)
        if Type>80
            error('dipcm:agarguments','The color index range is [1~80]');
        end
        Cmap = ColorsList{Type};
        colorName = ColorNames(Type);  % Color name
    else
        Cpos   = strcmpi(Type, ColorNames);
        Cmap = ColorsList{Cpos};
        colorName = Type;          % Color name
    end
    
    Color= interpCMap(Cmap,num);
else
    Color = ColorNames;           % Color Names List
    dispallBars(Color,1,80,'dip'); % show all bars
end
end