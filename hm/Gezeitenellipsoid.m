% -------------------------------------------------------------------------
% Gezeitenellipsoid.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Form des Gezeitenellipsoids.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Vergleich strenges Ellipsoid mit Legendre Polynom und Näherung

theta1  =  linspace(0,2*pi,500);
phi1    =  linspace(0,2*pi,500);
exz   =  -0.1;
r1    =  1 - 2*exz*legendreP(2,cos(theta1))/3;
r2    =  1 - exz*cos(2*theta1);
r3    =  0*r2+1;
figure()
% h(1) = polarplot(theta1,r1-r2,'Color',Colors(2,:),'LineWidth',2);
h(1) = polarplot(theta1,r1,'Color',Colors(2,:),'LineWidth',2);
hold on;
h(2) = polarplot(theta1,r2,'Color',Colors(3,:),'LineWidth',2);
h(3) = polarplot(theta1,r3);
legend(h,'Exakt','Näherung','Kreis','location','best','FontSize', 14);
legend box off

% Kugel und Ellipsoid
% [theta,phi] = ndgrid(linspace(0,180,180),linspace(0,360,180));
[theta,phi] = ndgrid(linspace(0,180,30),linspace(0,360,30));
R        = 1;
h        = 2*exz*2*exz*legendreP(2,cos(theta))/3;

% Kugel
xk  = R.*sind(theta).*cosd(phi);
yk  = R.*sind(theta).*sind(phi);
zk  = R.*cosd(theta);

% Ellipsoid
xe  = (R+h).*sind(theta).*cosd(phi);
ye  = (R+h).*sind(theta).*sind(phi);
ze  = (R+h).*cosd(theta);

% % 2D Darstellung Verformung Wasserhuelle durch Gezeiten
% figure('Name',' Verformung Wasserhuelle durch Gezeiten');
% hold on
% plot(ze,ye,'Color', Colors(3,:),'Linewidth',3);
% plot(R.*cosd(phi(1,:)),R.*sind(phi(1,:)),'Color', Colors(2,:),'Linewidth',2);
% grid on;
% axis equal
% ylim([-1.5 1.5]);
% xlim([-1.5 1.5]);
% xlabel('a.u.');
% ylabel('a.u.');
% set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

% 3D Darstellung Verformung Wasserhuelle durch Gezeiten
figure('Name',' Verformung Wasserhuelle durch Gezeiten');
mesh(ze,ye,xe,gradient(ze));
hold on
axis equal
set(gca, 'Fontsize', 14, 'linewidth', 1);


% figure('Name',' Verformung Wasserhuelle durch Gezeiten');
% % surf(xk,yk,zk)
% surf(ze,ye,xe,gradient(ze))

% Volumen Kugel
VK = 4*pi*R^3/3

[xs,ys,zs] = sphere(100);
xs = R*xs(:);
ys = R*ys(:);
zs = R*zs(:);
P = [xs,ys,zs];
P = unique(P,'rows');
shp = alphaShape(P,1.5);
VKnum = volume(shp)

% Volumen Ellipsoid
xes = xe(:);
yes = ye(:);
zes = ze(:);
P = [zes,xes,yes];
P = unique(P,'rows');
shp = alphaShape(P,1.5);
% % Plot an alpha shape using an alpha radius of 1.5.
% figure()
% plot(shp)
VEnum = volume(shp)

% Absolute Groessen

MM = 7.346e22; 
RE = 6.371e06;
ME = 5.972e24;
MS = 1.9884e30;
rEM = 384.4e06;
rES = 149.6e09;

% Gezeitenbeschleunigung Mond
gM = 1.5*MM*RE^3/ME/rEM^3

% Gezeitenbeschleunigung Sonne
gS = 1.5*MS*RE^3/ME/rES^3

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------