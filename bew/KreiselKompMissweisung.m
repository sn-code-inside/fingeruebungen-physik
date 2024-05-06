% -------------------------------------------------------------------------
% KreiselKompMissweisung.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Missweisung des Kreiselkompasses
% 
% Beispiel zur Missweisung des Kreiselkompasses bei einer Bewegung
% auf einem Laengenkreis
%
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", "--", "-", ":"];

%% Variablen: Traegheitsmomente, Erdrotation, Position auf der Erde, ...
%
L1 = -40; % L1'''' aus der Aufgabenstellung in kg*m²/s
J23  = 0.02;
OmegaE = 2.*pi/(24*3600); % Winkelgeschwindigkeit der Erdrotation
bg = 20;
theta = (90.-bg)/180*pi; % Breitengrad
R = 6378e3;

% Maximale betrachtete Geschwindigkeit
vmax = 300; % im m/s
% Anzahl der Werte
nmax = 200;

for k=1:nmax
    v(k) = (k-1)*vmax/(nmax-1);
    OmegaL = v(k)/R;
    BetaBeschlWerte = @(x) BetaBeschl(x,theta,OmegaL,OmegaE,J23,L1);
    betamin(k) = fzero(BetaBeschlWerte,-0.1);
end

%% Abbildung erstellen
%
figure();
subplot(1, 2, 1);
plot(v(:),rad2deg(betamin(:)),'Linewidth',2,'Color',Colors(1,:), ...
    'LineStyle',Style{1});
axis([0 20 -3 0])
xlabel('{\it v} in m/s','FontSize',14);
ylabel('\beta \rm in °','FontSize',14);
grid on;
set(gca,'FontSize',16);
hold off;
subplot(1, 2, 2);
plot(v(:),rad2deg(betamin(:)),'Linewidth',2,'Color',Colors(1,:), ...
    'LineStyle',Style{1});
%axis([0 1500 -40 0])
xlabel('{\it v} in m/s','FontSize',14);
ylabel('\beta \rm in °','FontSize',14);
grid on;
set(gca,'FontSize',16);
hold off;
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktion deren Nullstelle gesucht wird
function ddbeta = BetaBeschl(beta,theta,OmegaL,OmegaE,J23,L1)
     ddbeta = -J23*sin(theta)*OmegaL*OmegaE ...
     - L1*(sin(theta)*sin(beta)*OmegaE+cos(beta)*OmegaL) ...
     - J23*sin(theta)^2*sin(beta)*cos(beta)*OmegaE^2 ...
     - J23*sin(theta)*(cos(theta)^2-sin(theta)^2)*OmegaE*OmegaL ...
     + J23*cos(beta)*sin(beta)*OmegaL^2;
end 
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
