function [colorList,colorName] = slanCM(varargin)
% =================================================================
% type : type of colomap
% num : number of colors
% -------------------------------------------------------------------------
% From : matplotlab | https://matplotlib.org/
% + Perceptually Uniform Sequential£º
%    viridis  plasma  inferno  magma  cividis
% + Sequential£º
%    Greys  Purples  Blues  Greens  Oranges  Reds  YlOrBr  YlOrRd  
%    OrRd  PuRd  RdPu  BuPu  GnBu  PuBu  YlGnBu  PuBuGn  BuGn  YlGn
%    binary  gist_yarg  gist_gray  gray  bone  pink  spring  summer
%    autumn  winter  cool  Wistia  hot  afmhot  gist_heat  copper
% + Diverging£º
%    PiYG  PRGn  BrBG  PuOr  RdGy  RdBu  RdYlBu  RdYlGn  Spectral  coolwarm  bwr  seismic
% + Cyclic£º
%    twilight  twilight_shifted  hsv
% + Qualitative£º
%    Pastel1  Pastel2  Paired  Accent  Dark2  Set1  Set2  Set3  tab10  tab20  tab20b  tab20c
% + Miscellaneous£º
%    flag  prism  ocean  gist_earth  terrain  gist_stern  gnuplot  gnuplot2   CMRmap  
%    cubehelix  brg  gist_rainbow  rainbow  jet  turbo  nipy_spectral  gist_ncar
% -------------------------------------------------------------------------
% From : scicomap | https://github.com/ThomasBury/scicomap
% + diverging:
%    berlin  bjy  bky  BrBG  broc  bwr  coolwarm  curl  delta  fusion  guppy  iceburn  lisbon
%    PRGn  PiYG  pride  PuOr  RdBu  RdGy  RdYlBu  RdYlGn  redshift  roma  seismic  spectral
%    turbo  vanimo  vik  viola  waterlily  watermelon  wildfire
% + sequential:
%    afmhot  amber  amp  apple  autumn  batlow  bilbao  binary  Blues  bone  BuGn  BuPu
%    chroma  cividis  cool  copper  cosmic  deep  dense  dusk  eclipse  ember  fall  gem
%    gist_gray  gist_heat  gist_yarg  GnBu  Greens  gray  Greys  haline  hawaii  heat  hot
%    ice  inferno  imola  lapaz  magma  matter  neon  neutral  nuuk  ocean  OrRd  Oranges
%    pink  plasma  PuBu  PuBuGn  PuRd  Purples  rain  rainbow  rainbow-sc  rainforest  RdPu
%    Reds  savanna  sepia  speed  solar  spring  summer  tempo  thermal  thermal-2  tokyo
%    tropical  turbid  turku  viridis  winter  Wistia  YlGn  YlGnBu  YlOrBr  YlOrRd
% + circular:
%    bukavu  fes  infinity  infinity_s  oleron  rainbow-iso  romao  seasons  seasons_s
%    twilight  twilight_s
% + miscellaneous:
%    oxy  rainbow-kov  turbo
% + qualitative:
%    538  bold  brewer  colorblind  glasbey  glasbey_bw  glasbey_category10
%    glasbey_dark  glasbey_hv  glasbey_light  pastel  prism  vivid
% -------------------------------------------------------------------------
% From : cmasher | https://cmasher.readthedocs.io/
% + sequential:
%    amber  amethyst  apple  arctic  bubblegum  chroma  cosmic  dusk  eclipse  ember  emerald
%    fall  flamingo  freeze  gem  ghostlight  gothic  horizon  jungle  lavender  lilac  neon
%    neutral  nuclear  ocean  pepper  rainforest  sapphire  savanna  sepia  sunburst  swamp
%    torch  toxic  tree  tropical  voltage
% + diverging:
%   copper  emergency  fusion  guppy  holly  iceburn  infinity  pride  prinsenvlag  redshift
%   seasons  seaweed  viola  waterlily  watermelon  wildfire
% -------------------------------------------------------------------------
% https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap
% -------------------------------------------------------------------------
%   INPUT:
%      type:  Colormap type can be string or index number
%      num:  Color scale or level, default is 256
%
%   OUTPUT:
%      Color: A specified color matrix (double)
%      colorName: A specified color name (string)
%
%  Example:
%	   slanCM; %Show colorbar image of all color
%      [Color,colorName] = slanCM(60); 
% 	   [Color,colorName] = slanCM(60,256);
% 	   [Color,colorName] = slanCM(60,20); 
% --------------------------------------------------------
slanCM_Data = load('.\NIRSpecoctaneDemo\colorData\slanCM_Data.mat');
CList_Data = [slanCM_Data.slandarerCM(:).Colors];

narginchk(0,2)
if nargin<2
    num=-1;
else
    num = varargin{2};
end

if nargin<1
    colorName = slanCM_Data.fullNames;
    %disp(ctype);    
    dispallBars(colorName,121,200,'slan'); % Show all bars
    return
end

type = varargin{1};
if isnumeric(type)
    if type>200
        error('slancm:agarguments','The color index range is [1~200]');
    end
    Cmap = CList_Data{type};
    colorName =slanCM_Data.fullNames{type};
else
    Cpos  = strcmpi(type, slanCM_Data.fullNames);
    Cmap = CList_Data{Cpos};
    colorName = type;
end

% Interpolate to the setted color numbers
if num>0
    colorList = interpCMap(Cmap, num); 
else
    colorList = Cmap;    
end
end