% -------------------------------------------------------------------------
% Gradient.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet stellt Funktionen und ihren Gradienten dar.
%
% Nach: 
% Javier Hasbun (2024). 
% Classical_Mechanics_with_Matlab_applications_JEH.zip 
% https://www.mathworks.com/matlabcentral/fileexchange/20579-classical_mechanics_with_matlab_applications_jeh-zip), 
% MATLAB Central File Exchange. 
% Retrieved March 12, 2024.
%
%
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')

axmax=2.0;
xmax=axmax; ymax=axmax; zmax=axmax; % x,y,z Limits

vs=0.1;
xs=vs; ys=vs; zs=vs;                % Inkremente
N=2*axmax/vs;                       % Anzahl Plotpunkte NxNxN
dv=0.1;
dx=dv; dy=dv; dz=dv; 


%% Funktion 1
m=round(N/2+1); zm=-zmax+(m-1)*zs;         % z-Wert für Plot von Phi(x,y,z)
[x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);

Phi=10*x.*y.*exp(-(x.^2+y.^2+z.^2));       % Phi(x,y,z)
[dfx,dfy,dfz] = gradient(Phi,dx,dy,dz);    % Gradient von Phi(x,y,z)

mesh(x(:,:,m),y(:,:,m),Phi(:,:,m)) 

figure(1)
surf(x(:,:,m),y(:,:,m),Phi(:,:,m)) % 3D Plot
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
zlabel('Phi(x,y,z)','FontSize',14)
str=cat(2,'Phi(x,y,z)=10*x*y*exp(-(x^2+y^2+z^2)) @ ','z=',num2str(zm,3));
axis equal
h=title(str);
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',12);

figure(2)
[C,h] = contour(x(:,:,m),y(:,:,m),Phi(:,:,m)); % Contour Plot
clabel(C,h,'FontSize',12)                      % Contour Labels
hold on
quiver(x(:,:,m),y(:,:,m),dfx(:,:,m),dfy(:,:,m))% Gradient als Pfeile
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
h=title(str);
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',12);
axis equal


%% Funktion 2 Yukawa-Potential

% C     = 10;
% alpha = 0.5;
% Phi = C*exp(-alpha*sqrt(x.^2+y.^2+z.^2))./sqrt(x.^2+y.^2+z.^2); %Phi(x,y,z)
% [dfx,dfy,dfz] = gradient(Phi,dx,dy,dz);           % Gradient von Phi(x,y,z)
% 
% m=1; zm=-zmax+(m-1)*zs;         % z-Wert für Plot von Phi(x,y,z)
% [x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);
% mesh(x(:,:,m),y(:,:,m),Phi(:,:,m)) 
% 
% figure(3)
% surf(x(:,:,m),y(:,:,m),Phi(:,:,m)) % 3D Plot
% xlabel('x','FontSize',14)
% ylabel('y','FontSize',14)
% zlabel('Phi(x,y,z)','FontSize',14)
% str=cat(2,'Yukawa-Potential');
% h=title(str);
% set(gca,'FontSize',14);
% set(h,'FontWeight','normal','FontSize',12);
% 
% figure(4)
% [C,h] = contour(x(:,:,m),y(:,:,m),Phi(:,:,m)); % Contour Plot
% clabel(C,h,'FontSize',12)                      % Contour Labels
% hold on
% quiver(x(:,:,m),y(:,:,m),dfx(:,:,m),dfy(:,:,m))% Gradient als Pfeile
% xlabel('x','FontSize',14)
% ylabel('y','FontSize',14)
% h=title(str);
% set(gca,'FontSize',14);
% set(h,'FontWeight','normal','FontSize',12);
% axis equal
% 
% %% Funktion 3 Zentralkraft-Potential mit Drehimpuls (effektives Potential)
% 
% C     = 1;
% L     = 1.10;
% r = sqrt(x.^2+y.^2+z.^2);
% Phi = -C./r + L^2./r.^2 ;                       %Phi(x,y,z)
% [dfx,dfy,dfz] = gradient(Phi,dx,dy,dz);     % Gradient von Phi(x,y,z)
% 
% m=1; %zm=-zmax+(m-1)*zs;         % z-Wert für Plot von Phi(x,y,z)
% [x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);
% mesh(x(:,:,m),y(:,:,m),Phi(:,:,m)) 
% 
% figure(5)
% surf(x(:,:,m),y(:,:,m),Phi(:,:,m)) % 3D Plot
% xlabel('x','FontSize',14)
% ylabel('y','FontSize',14)
% zlabel('Phi(x,y,z)','FontSize',14)
% str=cat(2,'Effektives Kepler-Potential');
% h=title(str);
% set(gca,'FontSize',14);
% set(h,'FontWeight','normal','FontSize',12);
% 
% figure(6)
% [C,h] = contour(x(:,:,m),y(:,:,m),Phi(:,:,m)); % Contour Plot
% clabel(C,h,'FontSize',12)                      % Contour Labels
% hold on
% quiver(x(:,:,m),y(:,:,m),dfx(:,:,m),dfy(:,:,m))% Gradient als Pfeile
% xlabel('x','FontSize',14)
% ylabel('y','FontSize',14)
% h=title(str);
% set(gca,'FontSize',14);
% set(h,'FontWeight','normal','FontSize',12);
% axis equal
