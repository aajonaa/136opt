clc; clear; close all;
data = load('NIRspectra_data.mat');
NIRspec = data.NIR';                               % IR Spectrum(Transposed)
octaneSamples = data.octane';               % Octane samples
dataNums =  size(NIRspec,1);                 % Data length
N = size(octaneSamples,2);                     % Sample Numbers

lambda = linspace(900,1700,dataNums)'; % Wavelength
X =  repmat(octaneSamples,[dataNums 1]);    % x-axis/Octane
Y =  repmat(lambda,[1,N]);                     % y-axis/Wavelength
Z =  NIRspec(1:dataNums,:);                   % z-axis/Absorption

% Color scheme
% C = colormap(hsv(N));
% C = colormap(idipCM(53,N));          % 80 colors
C = colormap(slanCM(150,N));      % 200 colors
% C = colormap(nclCM(490,N));        % 487 colors
% C = colormap(othercolor(200,N));  % 404 colors

% Draw 2D curve
for k = 1:N
    plot(Y(:,k),Z(:,k)); hold on   
end
title(' NIR Absorption Spectra (60 samples)','FontSize',12);
grid on,box on
xlabel('Wavelength (nm)')
ylabel('Absorption(au)')

% Background Color
set(gcf,'Color', [1 1 1])

figure,
% Draw 3D curve
for k = 1:N
    plot3(X(:,k),Y(:,k),Z(:,k),'Color',C(k,:),...
        'LineWidth',1.2);hold on
   %fill3(X(:,k),Y(:,k),Z(:,k),C(k,:),'FaceAlpha',0.8,'EdgeColor','none');  
   stem3(X(:,k),Y(:,k),Z(:,k),'Marker','none','Color',C(k,:));
end
title(' NIR Absorption Spectra (60 samples)','FontSize',12);
grid on,box on
axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) min(Z(:)) max(Z(:))])
%set(gca,'YDir','reverse')
xlabel('Samples(Octane)','Rotation', -21)
ylabel('Wavelength(nm)','Rotation', 21)
zlabel('Absorption(au)')
view([45, 30])

% Background Color
set(gcf,'Color', [1 1 1])