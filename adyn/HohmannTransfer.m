% -------------------------------------------------------------------------
% HohmannTransfer.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Geschwindigkeitsschübe für den Hohmannn-Transfer
% als Funktion der Verhältnisse von Anfangs- und Endradius der kreisförmigen
% Erdumlaufbahn und die Transferzeiten für interplanetrae Hohman-Transfers
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Parameter
GME = 398600.4415;   % G*ME in km^3/s^2
RE  = 6378;          % Erdradius in km
AE  = 1.4967e8;      % Astronomische Einheit in km
h1  = 600;           % Höhe Bahn
h2  = 650;
r1  = h1+RE;         % geozentrischer Radius
r2  = h2+RE;


%% Berechnung Geschwindigkeitszuwachs

v1 = sqrt(GME/r1);
v2 = sqrt(GME/r2);
aH = 0.5*(r1+r2);         % Transferellipse Halbachse
vP = sqrt(GME*r2/aH/r1);  % v bei Perigäum
vA = sqrt(GME*r1/aH/r2);  % v bei Apogäum

% Geschwindigkeitszuwachs
Dv1 = vP-v1;
Dv2 = v2-vA;
Dv  = Dv1 + Dv2;
Dvc = v2-v1;


%% Berechnung Geschwindigkeitszuwachs für Graphik

r1 = 1*RE;
r2 = r1*linspace(1,10000,100000);
v1 = sqrt(GME/r1);
v2 = sqrt(GME./r2);
aH = 0.5*(r1+r2);
vP = sqrt(GME*r2./aH/r1);
vA = sqrt(GME*r1./aH./r2);
Dv1 = (vP/v1-1);
Dv2 = (v2-vA)/v1;
Dv  = (Dv1 + Dv2);
vesc = sqrt(2)-1;
figure()
semilogx(r2/r1,Dv, 'linewidth',2);
hold on
semilogx(r2/r1,Dv1, 'linewidth',2);
semilogx(r2/r1,Dv2, 'linewidth',2);
semilogx(r2/r1,vesc*r2(:)./r2(:), 'linewidth',2, 'Linestyle', ':');
grid on
ylabel('\Deltav / v_1')
xlabel('r_2 / r_1')
legend('\Deltav_g / v_1','\Deltav_1 / v_1','\Deltav_2 / v_1',...
       'v_{esc} / v_1','location','east');   
legend box off
title('Hohmann-Transfer \Deltav')
set(gca,'FontSize',14);


%% Berechnung Transferzeiten 

G   = 6.671e-11;         % G in in m^3/s^2/kg
MS  = 1.989e30;          % Masse Sonne in kg 
GMS = G*MS*1e-9;         % G*MS in km^3/s^2
r1 = AE;  %Parkbahn Erde;
r2 = r1*linspace(1,101,1001);
aH = (r1+r2)/2;
tH = pi*sqrt(aH.^3/GMS);
figure()
loglog(r2/AE,tH/86400/365, 'linewidth',2);
grid on
ylabel('t_H in Jahren')
xlabel('r_2 in AE')
title('Hohmann-Transferzeiten')
set(gca,'FontSize',14);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
