% -------------------------------------------------------------------------
% Rotation.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Vektorfunktion und ihre Rotation.
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
Colors = GetColorLines;

vmin=-1; vmax=1;
xmax=vmax; ymax=vmax; zmax=vmax;    % x,y,z Limits
xmin=vmin; ymin=vmin; zmin=vmin;
vs=0.25; xs=vs; ys=vs; zs=vs;       % Inkrement
N=(vmax-vmin)/vs;                   % Punktezahl NxNxN

m=round(N/2+1); zm=zmin+(m-1)*zs;   % z Wert für Plots

[x,y,z]=meshgrid(xmin:xs:xmax,ymin:ys:ymax,zmin:zs:zmax);

%% Funktion auswählen

% Example 1
fx=0.5*y.*z.^2;fy=-0.5*x.*z.^2;fz=0*x;str1='Vektorfunktion 0.5*(y*z^2,-x*z^2,0)';
% Example 2
% fx=y.*z.^2;fy=(-x.*z.^2);fz=0;str1='Vektorfunktion (y*z^2,-x*z^2,0)';
% Example 3
% fx=x.*y;fy=y.*z;fz=z.*x;str1='Vektorfunktion (xy,yz,zx)';
% Example 4
%   fx=-y;fy=x;fz=0.*z;str1='Vektorfunktion (-y,x,0)';


%%

figure ('name', 'Vektorfunktion')
quiver3(x,y,z,fx,fy,fz,1.5,'color',Colors(3,:)); 
str2=cat(2,str1,' @ z=',num2str(zm,3));
h= title(str2,'FontSize',14);
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
zlabel('z','FontSize',14)
axis equal
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',14);


%%

figure ('name', 'Vektorfunktion und Zirkulation auf einer Fläche')
[curlx,curly,curlz,cav] = curl(x,y,z,fx,fy,fz);%Roatation in rad/sec von f
x2=x(:,:,m); y2=y(:,:,m); fx2=fx(:,:,m); fy2=fy(:,:,m);
cav2 = curl(x2,y2,fx2,fy2);% omega_z = Zirkulation = Durchschnitt der Winkelgeschwindigkeit
surf(x2,y2,cav2)  % Wert von omega_z auf der z-Fläche z=zm
shading interp
hold on 
quiver(x2,y2,fx2,fy2);    %Vektorfunktion @ z=zm
str2=cat(2,'\omega_z-Wert und Vektorfunktion auf Fläche ','z=',num2str(zm,3));
str4=cat(2,'\omega_z, f_x, f_y, ');
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
zlabel(str4,'FontSize',14)
h= title(str2,'FontSize',14);
axis equal
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',14);
hold off

%%

figure ('name', 'Rotation der Vektorfunktion')
quiver3(x,y,z,curlx,curly,curlz,2); %draws curl's result as arrows in three dimensions
str3=cat(2,'Rotation der Vektorfunktion');
title (str3,'FontSize',14) 
xlabel('x','FontSize',14), ylabel('y','FontSize',14)
zlabel('z','FontSize',14)
h= title(str3,'FontSize',14);
axis equal
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',14);
