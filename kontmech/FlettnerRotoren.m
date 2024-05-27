% -------------------------------------------------------------------------
% FlettnerRotoren.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger 
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den AutorenCartarius
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Leistung eines Flettner-Rotors
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')

L = 30;        % H�he des Rotors
R = 5;         % Radisu des Rotors
omega = 3;     % Umdrehungszahl des Rotors in 1/s ungef�hr 180 /min
v = 20;        % Quergeschwindigkeit Wind in m/s
rhoL = 1.204;  % Dichte Luft kg/m^3

% Vorw�rtskraft in N
F = 2*pi*rhoL*L*R^2*v*omega;
fprintf(' Vorw�rtskraft in N : %6.2e',F);
% Schiffsvortriebsgeschwindigkeit in m/s
vs  = 6;
% Leistung in kW
P = F*vs;
fprintf('\n Leistung in kW     : %6.2e \n',P/1000);


