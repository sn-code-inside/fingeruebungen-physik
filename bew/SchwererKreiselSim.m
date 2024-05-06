% -------------------------------------------------------------------------
% SchwererKreiselSim.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Beispielrechnungen und Simulation zur Präzession des 
% schweren symmetrischen Kreisels.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", "--", "-", ":"];

%% Variablen: Traegheitsmomente, Masse
%
J3     = 1.0;
J12    = 2*J3/3;
omega3 = 10.;
m      = 1.;
g      = 9.81;
s      = 0.2;

% Integrationszeit in s
tmax    = 5.;
% Komponenten der Rotationsachse im koerperfesten System zu t=0
phi0   =  0.;
psi0   =  0.;
% Erhaltungsgroessen
pphi    = 0.5;
omega3  = 1.;

%% Bestimme zwei Werte fuer Energien mit und ohne Nutation
% 
% Suche das Minimum fuer die reine Praezessionsbewegung
theta0 = fminsearch(@(x)EffektivesPotential(x,pphi,omega3,J12, ...
    J3,m,g,s),1.);
fprintf('\n theta_0 = ');
theta0d = rad2deg(theta0);
fprintf(num2str(theta0d,'%4.1f')');
fprintf('°\n');
dtheta0=0;  %Anfangsgeschwindigkeit 0
Eges_0 = (EffektivesPotential(theta0,pphi,omega3,J12,J3,m,g,s)...
       +0.5*J12*dtheta0^2)/J3/omega3^2;

% Waehle eine Anfangsgeschwindigkeit und damit eine hoehere Energie fuer
% eine Bewegung mit Nutation
dtheta1 = 2.0;

%% Zeitentwicklung theta
%
% Loese Bewegungsgleichung fuer Fall ohne Nutation
tspan = [0,tmax];
Y0=[theta0;dtheta0;phi0;psi0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,pphi,omega3,J12,J3,m,g,s);

% ... fuer Fall mit Nutation
Y0=[theta0;dtheta1;phi0;psi0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t2,Y2] = ode45(@DGL,tspan,Y0,options,pphi,omega3,J12,J3,m,g,s);

%% Erstelle Abbildungen
%
% ohne Nutation
figure
subplot(1,2,1)
% theta
plot(t1,rad2deg(Y1(:,1)),'Linewidth',1,'Color',Colors(11,:), ...
    'LineStyle',Style{1});
ylabel('\phi,\vartheta,\psi','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
% phi
plot(t1,rad2deg(Y1(:,3)),'Linewidth',1,'Color',Colors(4,:), ...
    'LineStyle',Style{1});
% psi
plot(t1,rad2deg(Y1(:,4)),'Linewidth',1,'Color',Colors(3,:), ...
    'LineStyle',Style{1});
legend('\phi','\vartheta','\psi in °','FontSize',14,...
       'Location','northwest');
legend box off
set(gca,'FontSize',14);
hold off
grid on;
axis([0 tmax 0 400]);
% Darstellung der Figurenachse im dreidimensionalen Raum
subplot(1,2,2)
plot3(sin(Y1(:,4)).*sin(Y1(:,1)),cos(Y1(:,4)).*sin(Y1(:,1)), ...
    cos(Y1(:,1)),'Linewidth',1,'Color',Colors(5,:),'LineStyle',Style(1));
zlabel('\it z','FontSize',14);
ylabel('\it y','FontSize',14);
xlabel('\it x','FontSize',14);
grid on
hold on
for k = 1:2:length(Y1(:,1))
plot3([0; sin(Y1(k,4)).*sin(Y1(k,1))], [0; cos(Y1(k,4)).*sin(Y1(k,1))], ...
    [0; cos(Y1(k,1))],'Color',Colors(5,:),...
      'LineStyle',Style(1));
end
axis([-1 1 -1 1 -1 1]);
legend off
set(gca,'FontSize',14);

% mit Nutation
figure()
subplot(1,2,1)
% theta
plot(t2,rad2deg(Y2(:,1)),'Linewidth',1,'Color',Colors(11,:), ...
    'LineStyle',Style{1});
ylabel('\phi,\vartheta,\psi in °','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
% phi
plot(t2,rad2deg(Y2(:,3)),'Linewidth',1,'Color',Colors(4,:), ...
    'LineStyle',Style{1});
% psi
plot(t2,rad2deg(Y2(:,4)),'Linewidth',1,'Color',Colors(3,:), ...
    'LineStyle',Style{1});
legend('\phi','\vartheta','\psi','FontSize',14,...
       'Location','northwest');
legend box off
set(gca,'FontSize',14);
hold off
grid on;
axis([0 tmax 0 400]);
% Darstellung der Figurenachse im dreidimensionalen Raum
subplot(1,2,2)
plot3(sin(Y2(:,4)).*sin(Y2(:,1)),cos(Y2(:,4)).*sin(Y2(:,1)), ...
    cos(Y2(:,1)),'Linewidth',1,'Color',Colors(5,:),'LineStyle',Style(1));
zlabel('\it z','FontSize',14);
ylabel('\it y','FontSize',14);
xlabel('\it x','FontSize',14);
grid on
hold on
for k = 1:2:length(Y2(:,1))
plot3([0; sin(Y2(k,4)).*sin(Y2(k,1))], [0; cos(Y2(k,4)).*sin(Y2(k,1))], ...
    [0; cos(Y2(k,1))],'Color',Colors(5,:),...
      'LineStyle',Style(1));
end
axis([-1 1 -1 1 -1 1]);
legend off
set(gca,'FontSize',14);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

% Differentialgleichung
function dYdt = DGL(t,Y,pphi,omega3,J12,J3,m,g,s)
    % Y(1) = theta
    % Y(2) = d(theta)/dt
    % Y(3) = phi
    % Y(4) = psi
    dYdt = [Y(2)
            -(pphi-J3*omega3*cos(Y(1)))*J3*omega3/(J12^2*sin(Y(1)))+...
            (pphi-J3*omega3*cos(Y(1)))^2*cos(Y(1))/(J12^2*sin(Y(1))^3)+...
            m*g*s*sin(Y(1))/J12
            (pphi-J3*omega3*cos(Y(1)))/(J12*sin(Y(1))^2)
            omega3-(pphi-J3*omega3*cos(Y(1)))*cos(Y(1))/(J12*sin(Y(1))^2)];
end

function Pot = EffektivesPotential(theta,pphi,omega3,J12,J3,m,g,s)
     Pot = (pphi - J3*omega3*cos(theta)).^2./(2.*J12*sin(theta).^2) ...
           + 0.5*J3*omega3^2 + m*g*s*cos(theta);
end 
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
