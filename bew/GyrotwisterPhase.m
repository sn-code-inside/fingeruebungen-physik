% -------------------------------------------------------------------------
% GyrotwisterPhase.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Anregen des Gyrotwisters durch Bewegung der Hand. In dieser Programmversion
% wird eine Phasenverschiebung der Anregung der Hand zum Idealfall
% betrachtet.
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
pd     = 0.1*pi; % Phasendifferenz zwischen omega1 und Anregung


%% Numerische Berechnung
tv      = linspace(0.,tmax,1000);
Y0      = [omega1;omega3;phi1];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 

% Phasendiferenz 0
pd = 0;
[t1,Y1] = ode45(@DGL,tv,Y0,options,J12,J3,omegaHamp,pd);

% Phasendiferenz 0.2*pi
pd = 0.2*pi;
[t2,Y2] = ode45(@DGL,tv,Y0,options,J12,J3,omegaHamp,pd);

%Phasendifferenz 0.5*pi
pd = pi/2.;
[t3,Y3] = ode45(@DGL,tv,Y0,options,J12,J3,omegaHamp,pd);

%Phasendifferenz 0pi
pd = pi;
[t4,Y4] = ode45(@DGL,tv,Y0,options,J12,J3,omegaHamp,pd);

%% Graphische Ausgabe
figure()
subplot(1,2,1)
plot(t1,Y1(:,2)/pi,'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1});
axis([0 tmax 45 55]);
ylabel('\omega_3/\pi','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on;
plot(t2,Y2(:,2)/pi,'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1});
hold on;
plot(t3,Y3(:,2)/pi,'Linewidth',1,'Color',Colors(3,:),'LineStyle',Style{1});
hold on;
plot(t4,Y4(:,2)/pi,'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style{1});
grid on;
legend('\Delta\phi = 0','\Delta\phi = 0.2\pi','\Delta\phi = \pi/2',...
       '\Delta\phi = \pi','FontSize',14,'Location','southwest');
set(gca,'FontSize',16);  
hold off;

subplot(1,2,2)
plot(t1,Y1(:,1)/pi,'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1});
axis([0 2 -5 5]);
ylabel('\omega_1/\pi','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on;
plot(t2,Y2(:,1)/pi,'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1});
hold on;
plot(t3,Y3(:,1)/pi,'Linewidth',1,'Color',Colors(3,:),'LineStyle',Style{1});
hold on;
plot(t4,Y4(:,1)/pi,'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style{1});
grid on;
set(gca,'FontSize',16);  
hold off;
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%% Differentialgleichngen
function dYdt = DGL(t,Y,J12,J3,omegaHamp,pd)
    % Y(1) = omega1
    % Y(2) = omega3
    % Y(3) = phi1
    dYdt = [J3/J12*sin(Y(3))*sin(Y(3)+pd)*omegaHamp*Y(2)+ ...
            cos(Y(3))*sin(Y(3))*sin(Y(3)+pd)^2*omegaHamp
            J12/J3*Y(1)*omegaHamp*sin(Y(3))*sin(Y(3)+pd)
            Y(1)];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

