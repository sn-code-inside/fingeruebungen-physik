% -------------------------------------------------------------------------
% Billardspiel.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Streuwinkel beim Billardspiel
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

v      = 1;
R      = 1;
b      = linspace(0,2*R,180);
thetaS = 2*acosd(b/2/R);
theta1  = atand(sind(thetaS)./(1+cosd(thetaS)));
theta2  = 90 - thetaS/2;

u1 = v*1./(cosd(theta1)+sind(theta1)./tand(theta2));
u2 = u1.*sind(theta1)./sind(theta2);

figure
title('Streuwinkel');
plot(b,thetaS,'Color',Colors(4,:),'Linewidth',2,'LineStyle',Style(1));
hold on
plot(b,theta1,'Color',Colors(2,:),'Linewidth',2,'LineStyle',Style(1));
plot(b,theta2,'Color',Colors(3,:),'Linewidth',2,'LineStyle',Style(1));
ylabel('Streuwinkel \theta´,\theta_1,\theta_2 in ° ','FontSize',14);
xlabel('Stoßparameter \it{b} \rm in \it{R} ','FontSize',14);
legend('\theta´ Schwerpunktsystem','\theta_1',...,
       '\theta_2');
legend box off
axis([0 2 0 180])
grid on
set(gca,'Fontsize', 16);


figure
title('Geschwindigkeiten');
plot(b,sqrt(u1.^2+u2.^2),'Color',Colors(4,:),'Linewidth',2,'LineStyle',Style(3));
hold on
plot(b,(v^2-u1.^2)/v^2,'Color',Colors(4,:),'Linewidth',2,'LineStyle',Style(1));
plot(b,u1,'Color',Colors(2,:),'Linewidth',2,'LineStyle',Style(1));
plot(b,u2,'Color',Colors(3,:),'Linewidth',2,'LineStyle',Style(1));
ylabel('Geschwindigkeiten in m/s ','FontSize',14);
xlabel('Stoßparameter \it{b} \rm in \it{R} ','FontSize',14);
legend('Gesamtenergie','Energietransfer','u_1','u_2','location','south');
legend box off
axis([0 2 0 1])
grid on
set(gca,'Fontsize', 16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

