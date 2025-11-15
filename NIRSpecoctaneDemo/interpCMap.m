function cmap = interpCMap(cmapIn, num)
% INTERPCMMP Interpolate colormap according to 
%   the specified number of colors
%
%   INPUT:
%      cmapIn:  Input colormap 
%      num:  Color scale or level.
%
%   OUTPUT:
%      cmap: Output colormap interpolated
%---------------------------------------------
     xIn = 0:(size(cmapIn, 1)-1);
     xOut  = linspace(0, xIn(end), num);
     cmap = interp1(xIn, cmapIn, xOut);
end