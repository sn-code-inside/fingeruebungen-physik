% -------------------------------------------------------------------------
% GradientFunktion1.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Funktionen und ihre Gradienten.
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

syms x y z Phi ...
     gradPhi  'real'
 
Phi(x,y,z) = 10*x.*y.*exp(-(x.^2+y.^2+z.^2));       % Phi(x,y,z)
vecxyz     = [x,y,z];
gradPhi    = gradient(Phi,vecxyz);


%% Auswertung der symbolischen Berechnung

axmax=2.0;
xmax=axmax; ymax=axmax; zmax=axmax; % x,y,z Limits

vs=0.1;
xs=vs; ys=vs; zs=vs;                % Inkremente
N=2*axmax/vs;                       % Anzahl Plotpunkte NxNxN
dv=0.1;
dx=dv; dy=dv; dz=dv; 

figure('name','Gradient der Funktion Phi (numerisch)')
m=round(N/2)+5; zm=-zmax+(m-1)*zs;  % z-Wert für Plot von Phi(x,y,z)
[x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);
rexp = exp(-x.^2-y.^2-z.^2);
dfx = 10*y.*rexp-20*x.^2.*y.*rexp;
dfy = 10*x.*rexp-20*x.*y.^2.*rexp;
dfz = -20*x.*y.*z.*rexp;            % Gradient von Phi(x,y,z)
PhiFunc = 10*x.*y.*exp(-(x.^2+y.^2+z.^2));
[C,h] = contour(x(:,:,m),y(:,:,m),PhiFunc(:,:,m)); % Contour Plot
clabel(C,h,'FontSize',12)                      % Contour Labels
hold on
str=cat(2,'Phi(x,y,z)=10*x*y*exp(-(x^2+y^2+z^2)) @ ','z=',num2str(zm,3));
quiver(x(:,:,m),y(:,:,m),dfx(:,:,m),dfy(:,:,m))% Gradient als Pfeile
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
h=title(str);
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',12);
axis equal



%% Numerische Berechnungen für Funktion 1

m=round(N/2)+5; zm=-zmax+(m-1)*zs;         % z-Wert für Plot von Phi(x,y,z)
[x,y,z]=meshgrid(-xmax:xs:xmax,-ymax:ys:ymax,-zmax:zs:zmax);

Phi=10*x.*y.*exp(-(x.^2+y.^2+z.^2));       % Phi(x,y,z)
[dfx,dfy,dfz] = gradient(Phi,dx,dy,dz);    % Gradient von Phi(x,y,z)

% mesh(x(:,:,m),y(:,:,m),Phi(:,:,m)) 

figure('name','Funktion Phi(x,y,z)')
surf(x(:,:,m),y(:,:,m),Phi(:,:,m)) % 3D Plot
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
zlabel('Phi(x,y,z)','FontSize',14)
str=cat(2,'Phi(x,y,z)=10*x*y*exp(-(x^2+y^2+z^2)) @ ','z=',num2str(zm,3));
axis equal
h=title(str);
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',12);

figure('name','Gradient der Funktion Phi (numerisch)')
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

