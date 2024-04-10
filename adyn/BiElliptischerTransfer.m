% -------------------------------------------------------------------------
% BiElliptischerTransfer.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Geschwindigkeitsbedarfe 
% für einen bielleptischen Transfer.
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


%% Berechnung nach Geschwindigkeitsbedarfs

gamma = linspace(2,1000,999);  %Verhältnis r2/r1
r1    = 2*RE;
v1    = sqrt(GME/r1);

% Hohmann-Transfer
r2 = r1*gamma;
v2 = sqrt(GME./r2);
aH = 0.5*(r1+r2);
vP = sqrt(GME*r2./aH/r1);
vA = sqrt(GME*r1./aH./r2);
Dv1 = (vP-v1);
Dv2 = (v2-vA);
delv_H  = (Dv1 + Dv2);

% Bielliptischer Transfer
delv_B = v1*(sqrt(2)-1)*(1 + 1./sqrt(gamma));


% Graphik Geschwindigkeitsbedarf

figure(1)
semilogx(gamma,delv_H/v1, 'linewidth',2,'color',Colors(2,:));
hold on
semilogx(gamma,delv_B/v1, 'linewidth',2,'color',Colors(3,:));
grid on
ylabel('\Delta v / v_1')
xlabel('\gamma')
legend('Hohmann','Bielliptisch',...
        'location','northeast');   
legend box off
ttl=title('Bielliptisch vs. Hohmann ');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',16);


% ------------------------------------------------------------------------% 
%-------------------------------------------------------------------------%
% Ende Programm
% ------------------------------------------------------------------------%


