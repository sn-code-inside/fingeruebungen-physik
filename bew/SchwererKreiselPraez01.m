% -------------------------------------------------------------------------
% SchwererKreiselPraez01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Beispielrechnungen und Simulation zur Präzession des 
% schweren symmetrischen Kreisels.
% 
% 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = {"-";"--";":"};

% Variablen: Traegheitsmomente, Masse
J12 = 0.4;
J3  = 1.5*J12;
omega3 = 10.;
m   = 1.;
g   = 9.81;
s   = 0.2;

% Integrationszeit
tmax = 5.;
% Komponenten der Rotationsachse im koerperfesten System zu t=0
phi0 = 0.;
theta0 = pi/3.;
psi0 = 0.;
% Erhaltungsgroessen
pphi = 0.5;
omega3 = 1.;
dtheta0 = 1.5;

% Erstelle Graph des Potentials
[xx,yy]=fplot(@(x)EffektivesPotential(x,pphi,omega3,J12,J3,m,g,s)/...
              J3/omega3^2,[0.05 pi-0.05]);
figure
plot(xx,yy,'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1})
Eges = (EffektivesPotential(theta0,pphi,omega3,J12,J3,m,g,s)...
       +0.5*J12*dtheta0^2)/J3/omega3^2;
axis([0 pi 0 10.]);
ylabel('{\it V}_{\rm eff}(\vartheta)/{\it J}_3\omega_3^2','FontSize',14);
xlabel('\vartheta','FontSize',14);
hold on
% Zeichne Gesamtenergie ein
plot([0,pi],[Eges,Eges],'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1})
set(gca,'FontSize',14);
hold off

tspan=[0.,tmax];
Y0=[theta0;dtheta0;phi0;psi0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,pphi,omega3,J12,J3,m,g,s);

figure
plot(t1,Y1(:,1),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1});
set(gca,'FontSize',14);
ylabel('\vartheta','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold off



% Suche das Minimum fuer die reine Praezessionsbewegung und überprüfe es
theta0 = fminsearch(@(x)EffektivesPotential(x,pphi,omega3,J12,J3,m,g,s),1.);
disp(theta0);
dtheta0=0;
[xx,yy]=fplot(@(x)EffektivesPotential(x,pphi,omega3,J12,J3,m,g,s)/...
              J3/omega3^2,[0.05 pi-0.05]);
figure
plot(xx,yy,'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1})
Eges = (EffektivesPotential(theta0,pphi,omega3,J12,J3,m,g,s)...
       +0.5*J12*dtheta0^2)/J3/omega3^2;
axis([0 pi 0 10.]);
ylabel('{\it V}_{\rm eff}(\vartheta)/{\it J}_3\omega_3^2','FontSize',14);
xlabel('\vartheta','FontSize',14);
hold on
% Zeichne Gesamtenergie ein
plot([0,pi],[Eges,Eges],'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1})
set(gca,'FontSize',14);
hold off


tspan=[0.,tmax];
Y0=[theta0;dtheta0;phi0;psi0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,pphi,omega3,J12,J3,m,g,s);

figure
plot(t1,Y1(:,2),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1});
ylabel('\phi,\vartheta,\psi','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
plot(t1,Y1(:,3),'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1});
plot(t1,Y1(:,4),'Linewidth',1,'Color',Colors(3,:),'LineStyle',Style{1});
legend('\phi','\vartheta','\psi','FontSize',14,...
       'Location','northwest');
set(gca,'FontSize',14);
hold off
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


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
     Pot = (pphi - J3*omega3*cos(theta))^2/(2.*J12*sin(theta)^2) ...
           + 0.5*J3*omega3^2 + m*g*s*cos(theta);
end 
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
