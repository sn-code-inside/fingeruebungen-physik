% -------------------------------------------------------------------------
% KreiselKomp.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Einstellen der Nordrichtung beim gedaempften Kreiselkompass
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = {"-";"--";":"};

%% Variablen:
L1 = -40; % L1'''' aus der Aufgabenstellung in kg*m²/s
J23  = 0.02;
OmegaE = 2.*pi/(24*3600); % Winkelgeschwindigkeit der Erdrotation
bg = 50.9;
theta = (90.-bg)/180*pi; % Breitengrad
gamma = 0.001; % Daempfungskonstante (klein fuer z.B. Quecksilber)

% Integrationszeit
hmax = 4; % Zahl der Stunden, die betrachtet werden sollen
tmax = hmax*3600.;
% Startwert fuer beta, Start in Ruhe
beta1 = pi/6; % anfaengliche Misseisung
dbeta = 0.;


%% Numerische Berechnung 
tv      = linspace(0.,tmax,10000*hmax);
Y0      = [beta1;dbeta];
options = odeset('AbsTol',1.e-10,'RelTol',1.e-8); 
[t1,Y1] = ode45(@DGL,tv,Y0,options,L1,J23,OmegaE,theta,gamma);
beta2 = pi/10;
Y0      = [beta2;dbeta];
[t2,Y2] = ode45(@DGL,tv,Y0,options,L1,J23,OmegaE,theta,gamma);
beta3 = 2*pi/3;
Y0      = [beta3;dbeta];
[t3,Y3] = ode45(@DGL,tv,Y0,options,L1,J23,OmegaE,theta,gamma);


%% Graphische Ausgabe: Abbildung fuer eine Missweisung
figure()

subplot(1,2,1)
plot(t1/60.,Y1(:,1)/pi,'Linewidth',1,'Color',Colors(1,:),'LineStyle', ...
    Style{1});
axis([0 10 -beta1*1.05/pi beta1*1.05/pi]);
ylabel('\beta/\pi','FontSize',14);
xlabel('{\it t} \rm in min ','FontSize',14);
grid on;
set(gca,'FontSize',16);
hold off;

subplot(1,2,2)
plot(t1/3600.,Y1(:,1)/pi,'Linewidth',1,'Color',Colors(1,:),'LineStyle', ...
    Style{1});
axis([0 hmax -beta1*1.05/pi beta1*1.05/pi]);
ylabel('\beta/\pi','FontSize',14);
xlabel('{\it t} \rm in h ','FontSize',14);
grid on;
set(gca,'FontSize',16);
title('\beta(0) = \pi/6','FontWeight','normal');
hold off;

%% Abbildung zum vergleich verschiedener Missweisungen
figure()

subplot(1,2,1)
plot(t2/3600.,Y2(:,1)/pi,'Linewidth',1,'Color',Colors(1,:),'LineStyle', ...
    Style{1});
axis([0 hmax -beta2*1.05/pi beta2*1.05/pi]);
ylabel('\beta/\pi','FontSize',14);
xlabel('{\it t} \rm in h ','FontSize',14);
grid on;
set(gca,'FontSize',16);
title('\beta(0) = \pi/10','FontWeight','normal');
hold off;

subplot(1,2,2)
plot(t3/3600.,Y3(:,1)/pi,'Linewidth',1,'Color',Colors(1,:),'LineStyle', ...
    Style{1});
axis([0 hmax -beta3*1.05/pi beta3*1.05/pi]);
ylabel('\beta/\pi','FontSize',14);
xlabel('{\it t} \rm in h ','FontSize',14);
grid on;
set(gca,'FontSize',16);
title('\beta(0) = 2\pi/3','FontWeight','normal');
hold off;
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%% Differentialgleichng
function dYdt = DGL(t,Y,L1,J23,OmegaE,theta,gamma)
    % Y(1) = beta
    % Y(2) = d(beta)/dt
    dYdt = [Y(2)
            (L1/J23*sin(theta)*OmegaE+sin(theta)^2*OmegaE^2*cos(Y(1)) ...
            )*sin(Y(1))-gamma*Y(2)];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------