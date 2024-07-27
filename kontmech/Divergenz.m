% -------------------------------------------------------------------------
% Divergenz.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Vektorfunktionen und ihre Divergenz.
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

vmin=-1; vmax=1;
xmax=vmax; ymax=vmax; zmax=vmax; % x,y,z Limits
xmin=vmin; ymin=vmin; zmin=vmin;
vs=0.4; xs=vs; ys=vs; zs=vs; % Inkrement
N=(vmax-vmin)/vs; % Punktezahl NxNxN
m=round(N/2+1); zm=zmin+(m-1)*zs; % z Wert für Plotting
[x,y,z]=meshgrid(xmin:xs:xmax,ymin:ys:ymax,zmin:zs:zmax);
% Vektorfunktion ist f = fx i + fy j + fz k. 
% Verschiedene Funktionen zur Auswahl 
% Example 1
% fx=x-z; fy=x.*z; fz=y.^2;
% str1='Funktion1 vecf = (x-z,x*z,y^2)';
% Example 2
% fx=5*x.^2; fy=3*y.^2; fz=0.*z;
% str1='Funktion2 vecf = (5x^2,3y^2,0)';
% Example 3
fx=x; fy=y; fz=z; 
str1='Funktion3 vecf = (x,y,z)';
% Example 4
% fx=x.^2.*y; fy=y.^2.*z; fz=z.^2.*x;
% str1='Funktion4 vecf = (x.^2.*y,y.^2.*z,z.^2.*x))';

div = divergence(x,y,z,fx,fy,fz); %Divergenz der Vektorfunktion

figure ('name', 'Vektorfunktion')
quiver3(x,y,z,fx,fy,fz,1); 
str2=cat(2,str1,' @ z=',num2str(zm,3));
h= title(str2,'FontSize',14);
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
axis equal
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',14);


figure ('name', 'Divergenz')
surf(x(:,:,m),y(:,:,m),div(:,:,m)) 
str2=cat(2, str1, ' Divergenz @ z=',num2str(zm,3));
h= title(str2,'FontSize',14);
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
axis equal
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',14);
