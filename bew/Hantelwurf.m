% -------------------------------------------------------------------------
% Hantelwurf.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet ueber die numerische Loesung die Bahn eines Hantelwurfs 
% mit ODE45.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
P1.m1 = 0.1;
P1.m2 = 0.1;
P1.L =  1.0;
P1.g =  9.81;

P2.m1 = 0.1;
P2.m2 = 0.1;
P2.L =  1.0;
P2.g =  9.81;

v0  =  20;
alp =  65;
tend = 50;
v0x  = v0*cosd(alp);
v0y  = v0*sind(alp);
tspan = linspace(0,10,tend);
Y01 = [0; v0x; P1.L; v0y; -pi/2; 2];
Y02 = [0; v0x; P2.L; v0y; pi; 4];

% DGl. Solver
opts = odeset('Mass',@(t,Y) Mass(t,Y,P1));
[t,Y1] = ode45(@(t,Y) f(t,Y,P1),tspan,Y01,opts);
[t,Y2] = ode45(@(t,Y) f(t,Y,P2),tspan,Y02,opts);
  
xwMax = max(Y1(length(Y1(:,1)),1),Y2(length(Y1(:,1)),1));
xw1 = linspace(0, xwMax, 100);
zw1 = P1.L/2 + tand(alp)*xw1 - P1.g*xw1.^2/2/v0^2/(cosd(alp))^2;

zw2 = P2.L + tand(alp)*xw1 - P2.g*xw1.^2/2/v0^2/(cosd(alp))^2;
xw2 = xw1-P2.L/2;

figure
title('Trajektorie einer Hantel');
axis([-5 40 -0.05 20])
hold on
ylabel('Hoehe \it{z} \rm in m ','FontSize',14);
xlabel('Wurfweite \it{z} \rm in m ','FontSize',14);
grid on
xm1 = Y1(:,1) + 0.5*P1.L*cos(Y1(:,5));
zm1 = Y1(:,3) + 0.5*P1.L*sin(Y1(:,5));
plot(xm1,zm1,'Color',Colors(3,:),'Linewidth',1,'LineStyle',Style(2));
xm2 = Y2(:,1) + 0.5*P2.L*cos(Y2(:,5));
zm2 = Y2(:,3) + 0.5*P2.L*sin(Y2(:,5));
plot(xm2,zm2,'Color',Colors(9,:),'Linewidth',1,'LineStyle',Style(3));

plot(xw1,zw1,'Color',Colors(4,:),'Linewidth',1,'LineStyle',Style(2));
plot(xw2,zw2,'Color',Colors(4,:),'Linewidth',1,'LineStyle',Style(3));

for j = 1:length(t)
   theta1 = Y1(j,5);
   x1 = Y1(j,1);
   z1 = Y1(j,3);
   xvals = [x1 x1+P1.L*cos(theta1)];
   zvals = [z1 z1+P1.L*sin(theta1)];
   plot(xvals,zvals,'Color',Colors(2,:),'Linewidth',1);
   plot(xvals(1),zvals(1),'o','MarkerSize',7,...
    'MarkerEdgeColor',Colors(2,:),'MarkerFaceColor',Colors(2,:));
   plot(xvals(2),zvals(2),'o','MarkerSize',7,...
    'MarkerEdgeColor',Colors(3,:),'MarkerFaceColor',Colors(3,:));
   theta2 = Y2(j,5);
   x2 = Y2(j,1);
   z2 = Y2(j,3);
   xvals = [x2 x2+P2.L*cos(theta2)];
   zvals = [z2 z2+P2.L*sin(theta2)];
   plot(xvals,zvals,'Color',Colors(5,:),'Linewidth',1);
   plot(xvals(1),zvals(1),'o','MarkerSize',7,...
    'MarkerEdgeColor',Colors(5,:),'MarkerFaceColor',Colors(5,:));
   plot(xvals(2),zvals(2),'o','MarkerSize',7,...
    'MarkerEdgeColor',Colors(9,:),'MarkerFaceColor',Colors(9,:));
end
hold off
legend('Hantel von unten','Hantel ueber Kopf',...,
       'schraeger Wurf', 'schraeger Wurf','Location','south','NumColumns',2);
legend box off
set(gca,'Fontsize', 16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Funktionen 

function M = Mass(t,q,P)
% Parameter
m1 = P.m1;
m2 = P.m2;
L = P.L;
g = P.g;

% Massen Matrix
M = zeros(6,6);
M(1,1) = 1;
M(2,2) = m1 + m2;
M(2,6) = -m2*L*sin(q(5));
M(3,3) = 1;
M(4,4) = m1 + m2;
M(4,6) = m2*L*cos(q(5));
M(5,5) = 1;
M(6,2) = -L*sin(q(5));
M(6,4) = L*cos(q(5));
M(6,6) = L^2;
end

function dYdt = f(t,Y,P)
    % Parameter
    m1 = P.m1;
    m2 = P.m2;
    L = P.L;
    g = P.g;

    % DGl.
    dYdt = [Y(2)
        m2*L*Y(6)^2*cos(Y(5))
        Y(4)
        m2*L*Y(6)^2*sin(Y(5))-(m1+m2)*g
        Y(6)
        -g*L*cos(Y(5))];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
