% -------------------------------------------------------------------------
% NumericValettPendel.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Numerik Valett-Pendel
% -------------------------------------------------------------------------

%%
clc
clear 
close all 
addpath('IncludeFolder')
addpath('Data')
Colors=GetColorLines;

%% Initialisierung
% alle Daten in kg, m und s 

%%-------------------------------------------------------------------------
% Parameterfestlegung
% --------------------------------------------------------------------------

JT1 = 2.57e-3;
JT2 = 3.85e-3;
mF = 2;
DT = 0.018;
DT1 = DT*34./20.;
DT2 = DT*34./13.;
%DT2 = 
kappa1 = 6e-1;
kF = 11.6;

% Anfangswerte
Y10 = 0.0;              % phi1
Y20 = 0.0;              % dphi1  
Y30 = 0.0;              % phi2  
Y40 = 0.0;              % dphi2  
Y50 = 0.2;              % z
Y60 = 0;                % dz
Y0 = [Y10 Y20 Y30 Y40 Y50 Y60]; % Anfangsvektor

P1.JT1      = JT1;
P1.JT2      = JT2;
P1.mF       = mF;
P1.DT1      = DT1;
P1.DT2      = DT2;
P1.kappa1   = kappa1;
P1.kF       = kF;

%%--------------------------------------------------------------------------
% Lösen der Differentialgleichung
%--------------------------------------------------------------------------

% MATLABs Runge-Kutta ode45 Routine 
opts   = odeset('AbsTol',1.e-10,'RelTol',1.e-09);
ispan  = 501;    % Anzahl der berechneten Zeitschritte

tmax = 60;             
t = linspace(0,tmax,6000);

[t1,Y]      = ode45(@DGL_Valett_Lagrange,t, Y0,opts,P1);
phi1(:)     = Y(:,1);
phidot1(:)  = Y(:,2);
phi2(:)     = Y(:,3);
phidot2(:)  = Y(:,4);
z(:)        = Y(:,5);
zdot(:)     = Y(:,6);

% Energieberechnung
T(:) = 0.5*JT1.*phidot1(:).^2+0.5*JT2.*phidot2(:).^2+0.5*mF.*zdot(:).^2;
V(:) = 0.5*kF.*z(:).^2 + 0.5*DT1.*phi1(:).^2 + 0.5*DT2.*phi2(:).^2 ...
    - 0.5*kappa1.*z(:).*(phi1(:)-phi2(:));
En(:) = T(:) + V(:);

%%--------------------------------------------------------------------------
% Graphische Ausgabe
%--------------------------------------------------------------------------

figure(1)
subplot(311)
plot(t,z)
ylabel ('\it z \rm in m')
h=title('Lösung mit ODE45');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);


subplot(312)
plot(t,phi1/pi)
ylabel ('\phi_1 in \pi')

grid on
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(313),plot(t,phi2/pi)
xlabel('\it t \rm in s'), ylabel('\phi_2 in \pi'),
grid on
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

figure(2)
plot(t,En)
ylabel ('\it E \rm in J')
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

%%--------------------------------------------------------------------------
%  Funktionen
%--------------------------------------------------------------------------

function dY = DGL_Valett_Lagrange(t, Y, P1)
% solves lagrangian

dY     = [Y(2);
          -1/P1.JT1*(P1.DT1*Y(1)+P1.kappa1/2*Y(5));
          Y(4);
          -1/P1.JT2*(P1.DT2*Y(3)-P1.kappa1/2*Y(5));
          Y(6);
          -1/P1.mF*(P1.kappa1/2*Y(1)-P1.kappa1/2*Y(3)+P1.kF*Y(5))];
end