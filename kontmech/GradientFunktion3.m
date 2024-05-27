% -------------------------------------------------------------------------
% GradientFunktion3.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Funktionen und ihre Gradienten.
% Kepler-Potential mit Drehimpuls
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



%% Symbolische Berechnung

syms x y z r C L  alpha  ...
      'real'
 
PhiY(r)       = -C/r + L^2/r^2;                                %Phi(r)
Phixyz(x,y,z) = -C/sqrt(x^2+y^2+z^2) +L^2/(x^2+y^2+z^2);       % Phi(x,y,z)
vecxyz     = [x,y,z];
gradPhixyz = gradient(Phixyz,vecxyz)
gradPhir   = gradient(PhiY,r)


%% Auswertung der symbolischen Berechnung

axmax=2.0;
xmax=axmax; ymax=axmax; zmax=axmax; % x,y,z Limits

vs=0.1;
xs=vs; ys=vs; zs=vs;                % Inkremente
N=2*axmax/vs;                       % Anzahl Plotpunkte NxNxN
dv=0.1;
dx=dv; dy=dv; dz=dv; 

C     = 1;
L     = 1.10;

figure('name','Gradient der Funktion Phi (analytisch)')
m=round(N/2+19); zm=-zmax+(m-1)*zs;     % z-Wert für Plot von Phi(x,y,z)
[x1,y1,z1]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);

r1 = sqrt(x1.^2+y1.^2+z1.^2);
PhiFunc =  C./r1 - L^2./r1.^2;    

dfx = x1.*(C./r1.^3 - (2*L^2)./r1.^4);
dfy = y1.*(C./r1.^3 - (2*L^2)./r1.^4);
% Gradient von Phi(x,y,z)
[C,h] = contourf(x1(:,:,m),y1(:,:,m),PhiFunc(:,:,m)); % Contour Plot
clabel(C,h,'FontSize',12)                           % Contour Labels
hold on
str=cat(2,'Effektives Kepler-Potential');
quiver(x1(:,:,m),y1(:,:,m),dfx(:,:,m),dfy(:,:,m))% Gradient als Pfeile
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
h=title(str);
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',12);
axis equal

%% Funktion 3 Zentralkraft-Potential mit Drehimpuls (effektives Potential)

[x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);

[dfx,dfy,dfz] = gradient(PhiFunc,dx,dy,dz);% Gradient von Phi(x,y,z)

figure('name','Funktion Phi(x,y,z)')
% mesh(x(:,:,m),y(:,:,m),PhiFunc(:,:,m)) 

surf(x(:,:,m),y(:,:,m),PhiFunc(:,:,m)) % 3D Plot
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
zlabel('Phi(x,y,z)','FontSize',14)
h=title(str);
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',12);

figure('name','Gradient der Funktion Phi(x,y,z)')
[C,h] = contour(x(:,:,m),y(:,:,m),PhiFunc(:,:,m)); % Contour Plot
clabel(C,h,'FontSize',12)                      % Contour Labels
hold on
quiver(x(:,:,m),y(:,:,m),dfx(:,:,m),dfy(:,:,m))% Gradient als Pfeile
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
h=title(str);
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',12);
axis equal
