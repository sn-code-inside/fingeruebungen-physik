% -------------------------------------------------------------------------
% StokesFormel01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger 
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den AutorenCartarius
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet das Strömungsprofil für die Stokes-Reibung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

r0 = 1;  % Radius Kugel
v0 = 1;  % Geschwindigkeit im Unendlichen

z  = linspace(-4,4,81);
y  = linspace(-4,4,81);  
[z,y] = meshgrid(z,y);

%Berechnung Geschwindigkeitsfeld
theta = atan2(y,z);  % Winkel zwischen z-Richtung Geschwindkeit und 
                     % Positionsvektor

r = sqrt(z.^2+y.^2); % Abstand vom Zentrum der Kugel

for i =1:length(z)   % Kugeldimensionsausschluss
    for j=1:length(y)
             if r(i,j) < 0.99*r0 
             r(i,j) = NaN;
             end
    end
end

A = v0.*(1-(0.25*r0./r).*(3+1./r.^2));
B = -3.*v0*r0*(1-r0^2./r.^2)/4./r;

vz = A+B.*cos(theta).*cos(theta); %z-Komponente der Geschwindigkeit
vy = B.*cos(theta).*sin(theta);   %y-Komponente der Geschwindigkeit

figure ('name', 'Vektorfunktion')
ystart = [-0.2 -0.35 -0.5 -0.75 -1 -1.25 -1.5 -1.75 -2 -2.5 -3 -3.5 -4];
zstart = -4*ones(1,length(ystart)); 
h=streamline(z,y,vz,vy,zstart,ystart);
set(h,'color','red')
hold on
ystart=(-1)*[-0.2 -0.35 -0.5 -0.75 -1 -1.25 -1.5 -1.75 -2 -2.5 -3 -3.5 -4];
h=streamline(z,y,vz,vy,zstart,ystart);
%quiver(z,y,vz,vy,'color',Colors(3,:)); 
axis equal
xlabel('z','FontSize',14)
ylabel('r','FontSize',14)
str = "Stokes-Formel - Geschwindigkeitsfeld";
h = title(str);
axis equal
PlotCircle (0,0,1,Colors(2,:),4);
axis([-4 4 -4 4]);
grid on
set(gca,'FontSize',14);
set(h,'FontWeight','normal','FontSize',12);

%Berechnung Druck
eta = 1;
p = -3000/2*eta*v0*r0*z./r.^3;

figure()
contour(z,y,p,20,'ShowText','off','Linewidth',2); 
str = "Stokes-Formel - Druck";
h2 = title(str);
hold on
PlotCircle (0,0,1,Colors(2,:),4);
axis([-4 4 -4 4]);
xlabel('z','FontSize',14)
ylabel('r','FontSize',14)
axis equal

grid on
set(gca,'FontSize',14);
set(h2,'FontWeight','normal','FontSize',12);


