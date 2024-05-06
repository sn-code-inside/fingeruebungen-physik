% -------------------------------------------------------------------------
% .m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Dynamik des Nichtlinearen Pendels
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

Style = ["-", "-.", ":", "--", ":"];

L = 1;                          % Länge
g = 9.81;                       % Schwerebeschleunigung
omega0 = sqrt(g/L);             % Eigenfrequenz HO
tau0   = 2*pi/omega0;           % Schwingungsperiode HO
tmax   = 2*tau0;                % Zeitspanne
theta0 = 15;                    % Anfangsamplitude in °
thetaA  = deg2rad(theta0);

% Numerische Lösung
AB=[thetaA;0.0];       % AB für ode45
opt=odeset('AbsTol',1.e-8,'RelTol',1.e-5);      
[t,thetan]=ode45(@dgl_pendel,[0.0,tmax],AB,opt,omega0); % Numerische Lösung

% Ansatzmethode (Harmonische Näherung)
omega=omega0*sqrt(1-thetaA^2/8);
kappa=omega0^2/6;
A1=thetaA;
A3=kappa*A1^3/(27*kappa*A1^2-32*omega0^2);
theta1=thetaA*cos(omega*t)+A3*cos(3*omega*t);

% Störungsrechnung
omega  = omega0*(1-thetaA^2/16-7*thetaA^4/3072);
theta2 = thetaA*(1+thetaA^2/192+23*thetaA^4/(1024*36))*cos(omega*t) - ...
         thetaA^3/6*(1/32+thetaA^2/256)*cos(3*omega*t) + ...
         thetaA^5/1024/36*cos(5*omega*t);

% Harmonischer Oszillator
thetaHO=thetaA*cos(omega0*t);


%% Graphische Darstellung

figure()
plot(t,rad2deg(thetan(:,1)),'Color',Colors(2,:),'LineWidth',1,...
    'LineStyle', Style(1));
hold on;
plot(t,rad2deg(theta1),'Color',Colors(3,:),'LineWidth',1,'LineStyle',Style(2));
plot(t,rad2deg(theta2),'Color',Colors(4,:),'LineWidth',1,'LineStyle',Style(3));
plot(t,rad2deg(thetaHO),'Color',Colors(1,:),'LineWidth',1,'LineStyle',Style(1));
h=legend('Numerische Lösung','Harmon. Ansatz','Störungsrechnung',...
         'Harmon. Oszillator','NumColumns',2);
set(h,'FontSize',12)
axis([0 max(t) -rad2deg(thetaA) rad2deg(thetaA*(1+0.4))]);
str=cat(2,'\theta_A = ',num2str(rad2deg(thetaA),3),'°');
text(0.5,rad2deg(thetaA*(1+0.2)),str,'FontSize',12);
xlabel('Zeit in s','FontSize',14)
ylabel('Amplitude (°)','FontSize',14)
grid on;
legend box off;

figure()
plot(t,rad2deg(thetan(:,1)),'Color',Colors(2,:),'LineWidth',1,...
    'LineStyle', Style(1));
hold on;
plot(t,rad2deg(thetaHO),'Color',Colors(1,:),'LineWidth',1,'LineStyle',Style(1));
plot(t,10*(rad2deg(thetan(:,1))-rad2deg(thetaHO)),'Color',Colors(2,:),...
    'LineWidth',1,'LineStyle', Style(1));
h=legend('Numerische Lösung',...
         'Harmon. Oszillator','Abweichung*10','NumColumns',2);
set(h,'FontSize',12)
axis([0 max(t) -rad2deg(thetaA) rad2deg(thetaA*(1+0.4))]);
str=cat(2,'\theta_A = ',num2str(rad2deg(thetaA),3),'°');
text(0.5,rad2deg(thetaA*(1+0.2)),str,'FontSize',12);
xlabel('Zeit in s','FontSize',14)
ylabel('Amplitude (°)','FontSize',14)
grid on;
legend box off;
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% Differentialgleichung

function dY = dgl_pendel(t, Y, omega0)
    % Y(1)-Winkel
    % Y(2)-Winkelgeschwindigkeit
    dY = [Y(2); -omega0^2*sin(Y(1))];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
