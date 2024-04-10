% -------------------------------------------------------------------------
% PreleDraht.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Perle auf Draht
% 
% Programm berechnet die zeitliche Entwicklung der Perle auf Draht fuer
% verschiedene Parameter (ohne Reibung). 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

%%
% Anfangsbedingungen, Parameter
L=5;                                % Laenge Draht
m1=1;                               % Massen
g=9.81;                             % Schwerebeschleunigung
TP = 20;                            % Umlaufzeit 

% Kreisfrequenz Scheibe
omega0 = 2*pi/TP;  

tmax=2.5*TP;                        % Gesamtzeit und Schrittweite
steps = 100; 
tsim =(0.0:tmax/steps:tmax);
N=length(tsim);                     % Berechnungszeitraum

% Anfangsbedingungen Koordinate r
r0= 1; r_dot0 = 0.1;               
muS = 1.0;                          % StokesReibung Koeefizient in 1/s
phi0=deg2rad(45);                   % Anfangswinkel in rad
kf = 0.5*omega0^2;
kf = 0.2;

%%
% Berechnung Darstellung mit Stokes Reibung
% Loesung der Bewegungsgleichungen 

if kf < omega0^2
    lambda1   = - muS/2 + 0.5*sqrt(muS^2 - 4*(kf-omega0^2));
    lambda2   = - muS/2 - 0.5*sqrt(muS^2 - 4*(kf-omega0^2));
    A = (r_dot0-r0*lambda2)/(lambda1-lambda2);
    B = r0-A;
    rho    =  A*exp(lambda1*tsim) + B*exp(lambda2*tsim);
    rhodot =  A*lambda1*exp(lambda1*tsim) + B*lambda2*exp(lambda2*tsim);
    phi =  phi0 + omega0*tsim;          % in °

    rhodot(1) = r_dot0;
    rho(1) = r0;
    % Vorzeichenumkehr bei Nulldurchgang
    for k=1:length(rho)
        if rho(k) <0
            rho(k) = -rho(k);
            rhodot(k) = -rhodot(k);
            phi(k) = phi(k) + pi;
        end
    end
else
    lambda1   = - muS/2 + 0.5*sqrt(muS^2 - 4*(kf-omega0^2));
    lambda2   = - muS/2 - 0.5*sqrt(muS^2 - 4*(kf-omega0^2));
    A = (r_dot0-r0*lambda2)/(lambda1-lambda2);
    B = r0-A;
    rho =  exp(-muS.*tsim).*exp(lambda1.*tsim) + B*exp(lambda2.*tsim);
    phi =  phi0 + omega0*tsim;          % in °
    rhodot =  A*lambda1*exp(lambda1*tsim) + B*lambda2*exp(lambda2*tsim);
end


figure()
subplot(1,2,1)
semilogy(tsim,rho);
xlabel('Zeit in s','FontSize',12);
ylabel('Radius in m','FontSize',12);
str1 = 'Radius';
title(str1 ,'FontSize',12);
grid on

subplot(1,2,2)
semilogy(tsim,abs(rhodot));
xlabel('Zeit in s','FontSize',12);
ylabel('Radialgeschw. in m/s','FontSize',12);
grid on
str1 = 'Radialgeschwindigkeit';
title(str1 ,'FontSize',12);

%%
% Simulation

figure;
% Aufbereitung für Darstellung     
rho    =  rho/min(rho); % Normierung wegen Logarithmischer Darstellung
rmax=rho(N); 

str1 ='Perle auf Draht - Radius Logarithmisch';
p = polarplot([phi(1); phi(1)], [0; log(rmax)]);
p.Color = Colors(3,:);
p.LineStyle = '-.';
p.LineWidth = 1;
p.Visible = 'on';

rho = abs(rho);

hold on
% pbaspect([(ap-am)/1.2/L 1 1])           % Aspektverhältnis Achsen
for k=1:N
   % clf
   hold on
   % Masse m
   p(1) = polarplot(phi(k),log(rho(k)));
   p(1).Color = Colors(2,:);
   p(1).Marker = 'o';
   p(1).MarkerSize = 5;
   p(1).LineWidth = 2;
   p(1).Visible = 'on';
   p(2) = polarplot([phi(k); phi(k)], [0.001; log(rmax)]);
   p(2).Color = Colors(3,:);
   p(2).LineWidth = 1;
   p(2).Visible = 'on';
   pause(0.025)
   p(1).Visible = 'off';
   p(2).Visible = 'off';
end
p(2) = polarplot([phi(N); phi(N)], [0; log(rmax)]);
p(2).Color = Colors(3,:);
p(2).LineWidth = 1;
p(2).Visible = 'on';
p(1) = polarplot(phi(N),log(rho(N)));
p(1).Color = Colors(2,:);
p(1).Marker = '*';
p(1).MarkerSize = 5;
p(1).LineWidth = 2;
p(1).Visible = 'on';
hold on;
pt = polarplot(phi,log(rho),'Linewidth',1,'LineStyle','-.','color',Colors(3,:));
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------