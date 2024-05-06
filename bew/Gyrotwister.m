% -------------------------------------------------------------------------
% Gyrotwister.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Anregen des Gyrotwisters durch Bewegung der Hand.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = {"-";"--";":"};

%% Variablen:
J12  = 1e-5;
J3   = 2e-5;
omegaHamp = 0.05; % Bewegung der Hand, Amplitude

% Integrationszeit
tmax = 10.; % Zahl der Sekunden, die betrachtet werden sollen
% Startwert fuer die Winkelgeschwindigkeite und Winkel
omega1 = 0.5*2*pi; % 1/20 von omega3 (Abrollen mit Radienverhältnis 1:50)
omega3 = 25*2*pi; % Schnelle Rotation des Kreisels
phi1   = 0.; % Bewegung der Hand parallel zur Figurenachse zu Beginn


%% Numerische Berechnung
tv      = linspace(0.,tmax,1000);
Y0      = [omega1;omega3;phi1];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tv,Y0,options,J12,J3,omegaHamp);

%% Graphische Ausgabe
figure()
plot(t1,Y1(:,2)/pi,'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1});
axis([0 tmax 45 60]);
ylabel('\omega_3/\pi','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
grid on;
set(gca,'FontSize',16);  
hold off;
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Differentialgleichungen
function dYdt = DGL(t,Y,J12,J3,omegaHamp)
    % Y(1) = omega1
    % Y(2) = omega3
    % Y(3) = phi1
    dYdt = [J3/J12*sin(Y(3))^2*omegaHamp*Y(2)+cos(Y(3))*sin(Y(3))^3 ...
             *omegaHamp
            J12/J3*Y(1)*omegaHamp*sin(Y(3))^2
            Y(1)];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

