% -------------------------------------------------------------------------
% Kontraktion.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Kontraktion von Galaxien
% Energie einer rotierenden Galaxie
% 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
c    = 3e8;                     % Lichtgeschwindigkeit in m/s
MS   = 1.989e30;                % Sonnenmasse in kg
M    = 1e12*MS;                 % Masse Galaxie in Sonnemasse
Rend = 1e5;                     % Endradius in LJ (Androemda)
G    = 6.672e-11;               % Gravitationskonstante
lyf  = c*86400;
rend = Rend*lyf;                % Endradius in m

L    = sqrt(6*G*M^3*rend/25)    % Drehimpuls

rmin =rend/4;
q = linspace(rmin,10*rmin,1000); 

T= 5*L^2/4/M./q.^2;
U= -3*G*M^2/5./q;
qlj = q/lyf;

ymin = min(U);

%%
% Graphik
%
figure();
plot(qlj,U,'Color',Colors(3,:),'LineWidth',2,'LineStyle',Style(1));
hold on
plot(qlj,T,'Color',Colors(2,:),'LineWidth',2,'LineStyle',Style(1));
plot(qlj,T+U,'Color',Colors(4,:),'LineWidth',2,'LineStyle',Style(1));
% line([0,xmax],[0,0],'Color','k','LineStyle',Style(4));
axis([0, 10*rmin/lyf, 1.2*ymin, -1.2*ymin])
ylabel('Energie in J','FontSize',14);
xlabel('Radius in LJ','FontSize',14);
h=legend('Gravitation','Rotation','Gesamt'); set(h,'FontSize',14);
legend box off;
grid on
set(gca,'Fontsize', 16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------