% -------------------------------------------------------------------------
% Rutherford.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Wirkungsquerschnitte für
% Rutherford-Streuung mit und ohne Rückstoß
%
% Streuteilchen: alpha-Teilchen und Gold 
% Streuzentrum:  Gold
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
Marker = ['o','d','o','s','+'];

% Variablen
u     = 1.660539e-27;       % Atomare Masseneinheit kg
eq    = 1.602e-19;          % Elementarladung in As
eps0  = 8.854e-12;          % Dieelektrizitaetskonstante As/Vm

% Winkel im CoM-System
theta1deg  = linspace(30,88,30);
theta1     = deg2rad(theta1deg);

%% Gold auf Gold

ZT    = 79;                 % Ladung Streuteilchen
ZZ    = 79;                 % Ladung Streuzentrum
kappa  = eq^2*ZZ*ZT/4/pi/eps0;
mZ  =  197;                 % Masse Streuzentrum in u
mT  =  197;                 % Masse Streuteilchen in u
E      = 10e6;              % Energie im eV
E      = E/6.242e+18;       % Umrechnung in Ws
E_Lab  = E/((mT+mZ)/mZ);

% Winkel im Laborsystem
hdeg    = 0.05;             % step size
psideg = 30:hdeg:88;        % domain
psi    = deg2rad(psideg);
h      = deg2rad(hdeg);
gamma  = mT/mZ;
q      = 1./tan(psi);
ctheta = (1./(1+q.^2)).*(-gamma + sqrt(q.^2-gamma^2*q.^2+q.^4));
theta  = acos(ctheta);
dtheta = diff(theta)/h;
dtheta =[dtheta,dtheta(end)];

% Wirkungsquerschnitt korrekt berechnet allgemeine Verhältnisse
sigma_Au_Au = dtheta.*sin(theta)./...
              sin(psi).*(kappa/4/E_Lab)^2./sin(theta/2).^4;
% Wirkungsquerschnitt CoM, 
sigma_CoM_Au = (kappa/4/E)^2./sin(theta1/2).^4;
% Wirkungsquerschnitt Lab System gleiche Massen (mit Rückstoß)
sigma_EqM = 4.*cos(theta1).*(kappa/4/E_Lab)^2./sin(theta1).^4;


%% Silber auf Gold

mZ  =  197;                 % Masse Streuzentrum in u
mT  =  109;                 % Masse Streuteilchen in u
E      = 1e6;               % Energie im eV
E      = E/6.242e+18;       % Umrechnung in Ws
E_Lab  = E/((mT+mZ)/mZ);
ZT    = 4;                  % Ladung Streuteilchen
ZZ    = 97;                 % Ladung Streuzentrum
kappa  = eq^2*ZZ*ZT/4/pi/eps0;

% Winkel im Laborsystem
hdeg    = 0.05;             % step size
psideg = 30:hdeg:88;        % domain
psi    = deg2rad(psideg);
h      = deg2rad(hdeg);
gamma  = mT/mZ;
q      = 1./tan(psi);
ctheta = (1./(1+q.^2)).*(-gamma + sqrt(q.^2-gamma^2*q.^2+q.^4));
theta  = acos(ctheta);
dtheta = diff(theta)/h;
dtheta =[dtheta,dtheta(end)];

%Wirkungsquerschnitt CoM, 
sigma_CoM_Ag= (kappa/4/E)^2./sin(theta1/2).^4;
%Wirkungsquerschnitt korrekt berechnet allgemeine Verhältnisse
sigma_Ag_Au = dtheta.*sin(theta)./...
              sin(psi).*(kappa/4/E_Lab)^2./sin(theta/2).^4;

%Wirkungsquerschnitt CoM, 
sigma_CoM_Ag= (kappa/4/E)^2./sin(theta1/2).^4;
%Wirkungsquerschnitt korrekt berechnet allgemeine Verhältnisse
sigma_Ag_Au = dtheta.*sin(theta)./...
              sin(psi).*(kappa/4/E_Lab)^2./sin(theta/2).^4;

          
%% Graphische Ausgabe

figure();
semilogy(rad2deg(theta1),sigma_CoM_Au,'Color',Colors(2,:),'Linewidth',2);
hold on
semilogy(rad2deg(psi),sigma_Au_Au,'Color',Colors(3,:), 'Linewidth',2);
semilogy(rad2deg(theta1),sigma_CoM_Ag ,'Color',Colors(2,:),'Linewidth',2);
semilogy(rad2deg(psi),sigma_Ag_Au ,'Color',Colors(4,:),'Linewidth',2);
grid on
ax = gca;
ax.XLim = [30 90];
% ax.YLim = [1e-49 1e-47];
legend('CoM Au-Au','Laborsystem Au-Au','CoM Ag-Au','Laborsystem Ag-Au');
%      'Laborsystem \it m_T = m_Z');
legend box off
ylabel('d\sigma/d\Omega bzw. d\sigma/d\Omega^\prime in m^2','FontSize',14);
xlabel('\theta bzw. \theta^\prime','FontSize',14);
title('Au- und Ag-Teilchen an Goldfolie','FontWeight','normal');
set(gca,'FontSize',16)   


%% alpha auf Gold

mZ  =  197;                 % Masse Streuzentrum in u
mT  =  4;                   % Masse Streuteilchen in u
E      = 5e5;               % Energie im eV
E      = E/6.242e+18;       % Umrechnung in Ws
E_Lab  = E/((mT+mZ)/mZ);
ZT    = 2;                  % Ladung Streuteilchen
ZZ    = 97;                 % Ladung Streuzentrum
kappa  = eq^2*ZZ*ZT/4/pi/eps0;

% Winkel im Laborsystem
hdeg    = 0.05;             % step size
psideg = 30:hdeg:88;        % domain
psi    = deg2rad(psideg);
h      = deg2rad(hdeg);
gamma  = mT/mZ;
q      = 1./tan(psi);
ctheta = (1./(1+q.^2)).*(-gamma + sqrt(q.^2-gamma^2*q.^2+q.^4));
theta  = acos(ctheta);
dtheta = diff(theta)/h;
dtheta =[dtheta,dtheta(end)];

% Wirkungsquerschnitt CoM, 
sigma_CoM_Alpha_Au = (kappa/4/E)^2./sin(theta1/2).^4;
% Wirkungsquerschnitt korrekt berechnet allgemeine Verhältnisse
sigma_Alpha_Au  = dtheta.*sin(theta)./...
                  sin(psi).*(kappa/4/E_Lab)^2./sin(theta/2).^4;

%% alpha auf Silber

mZ  =  109;                 % Masse Streuzentrum in u
mT  =  4;                   % Masse Streuteilchen in u
E      = 5e5;               % Energie im eV
E      = E/6.242e+18;       % Umrechnung in Ws
E_Lab  = E/((mT+mZ)/mZ);
ZT    = 2;                  % Ladung Streuteilchen
ZZ    = 47;                 % Ladung Streuzentrum
kappa  = eq^2*ZZ*ZT/4/pi/eps0;

% Winkel im Laborsystem
hdeg    = 0.05;             % step size
psideg = 30:hdeg:88;        % domain
psi    = deg2rad(psideg);
h      = deg2rad(hdeg);
gamma  = mT/mZ;
q      = 1./tan(psi);
ctheta = (1./(1+q.^2)).*(-gamma + sqrt(q.^2-gamma^2*q.^2+q.^4));
theta  = acos(ctheta);
dtheta = diff(theta)/h;
dtheta =[dtheta,dtheta(end)];

% Wirkungsquerschnitt CoM, 
sigma_CoM_Alpha_Ag = (kappa/4/E)^2./sin(theta1/2).^4;
% Wirkungsquerschnitt korrekt berechnet allgemeine Verhältnisse
sigma_Alpha_Ag = dtheta.*sin(theta)./...
                 sin(psi).*(kappa/4/E_Lab)^2./sin(theta/2).^4;



%% Graphische Ausgabe

figure();
semilogy(rad2deg(theta1), sigma_CoM_Alpha_Au,'+','MarkerEdgeColor',...
         Colors(2,:),'MarkerSize',7, 'LineWidth',2);
hold on
semilogy(rad2deg(psi), sigma_Alpha_Au,'Color',Colors(3,:),'Linewidth',2);
semilogy(rad2deg(theta1), sigma_CoM_Alpha_Ag,'+','MarkerEdgeColor',...
         Colors(2,:),'MarkerSize',7, 'LineWidth',2);
semilogy(rad2deg(psi), sigma_Alpha_Ag,'Color',Colors(4,:),'Linewidth',2);
grid on
ax = gca;
ax.XLim = [30 90];
legend('CoM \alpha-Au','Laborsystem \alpha-Au','CoM \alpha-Ag',...
       'Laborsystem \alpha-Ag');
legend box off
ylabel('d\sigma/d\Omega bzw. d\sigma/d\Omega^\prime in m^2','FontSize',14);
xlabel('\theta bzw. \theta^\prime','FontSize',14);
title('\alpha-Teilchen an Au- und Ag-Folie','FontWeight','normal');
set(gca,'FontSize',16)   
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



