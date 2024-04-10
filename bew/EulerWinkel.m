% -------------------------------------------------------------------------
% EulerWinkel.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Euler-Winkel
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%%
OR = [0;0;0];


tk = linspace(0,2*pi,360);
xk(:,1)=cos(tk);
xk(:,2)=sin(tk);
xk(:,3)=0*tk;
t = linspace(0,1,10);
x1(:,1)=t;
x1(:,2)=0*t;
x1(:,3)=0*t;
y1(:,1)=0*t;
y1(:,2)=t;
y1(:,3)=0*t;
z1(:,1)=0*t;
z1(:,2)=0*t;
z1(:,3)=t;
x1 = x1';
y1 = y1';
z1 = z1';

phi   = deg2rad(-40);
theta = deg2rad(-35);
psi   = deg2rad(-15);
R1 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
R3 = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
R2 = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];

x2 = mtimes(R1,x1);
y2 = mtimes(R1,y1);
z2 = mtimes(R1,z1);
x3 = mtimes(R2,x2);
y3 = mtimes(R2,y2);
z3 = mtimes(R2,z2);
x4 = mtimes(R3,x3);
y4 = mtimes(R3,y3);
z4 = mtimes(R3,z3);

xk = xk';
xkf = mtimes(R1,xk);
xkf = mtimes(R2,xkf);
xkf = mtimes(R3,xkf);

figure()
line(x1(1,:),x1(2,:),x1(3,:),'Color',Colors(3,:),'LineWidth',2);
hold on
line(y1(1,:),y1(2,:),y1(3,:),'Color',Colors(3,:),'LineWidth',2);
line(z1(1,:),z1(2,:),z1(3,:),'Color',Colors(3,:),'LineWidth',2);
line(x2(1,:),x2(2,:),x2(3,:),'Color',Colors(2,:),'LineWidth',2);
line(y2(1,:),y2(2,:),y2(3,:),'Color',Colors(2,:),'LineWidth',2);
line(z2(1,:),z2(2,:),z2(3,:),'Color',Colors(2,:),'LineWidth',2);
line(x3(1,:),x3(2,:),x3(3,:),'Color',Colors(1,:),'LineWidth',2);
line(y3(1,:),y3(2,:),y3(3,:),'Color',Colors(1,:),'LineWidth',2);
line(z3(1,:),z3(2,:),z3(3,:),'Color',Colors(1,:),'LineWidth',2);
line(x4(1,:),x4(2,:),x4(3,:),'Color',Colors(4,:),'LineWidth',2);
line(y4(1,:),y4(2,:),y4(3,:),'Color',Colors(4,:),'LineWidth',2);
line(z4(1,:),z4(2,:),z4(3,:),'Color',Colors(4,:),'LineWidth',2);
line(xk(1,:),xk(2,:),xk(3,:),'Color',Colors(3,:),'LineWidth',2);
line(xkf(1,:),xkf(2,:),xkf(3,:),'Color',Colors(4,:),'LineWidth',2);

view(3)
xlabel('x in AE')
ylabel('y in AE');
zlabel('z in AE');
xlim([-2 2]);
ylim([-2 2]);
zlim([-2 2]);
xl = xlim;
yl = ylim;
[X,YB] = meshgrid(xl,yl);
surf(X,YB,zeros(size(X)))
shading flat
alpha 0.1

grid on
grid minor
set(gca,'FontSize',18);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

