% -------------------------------------------------------------------------
% TrebuchetLGL01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Floating Arm Trebuchet
% 
% Programm berechnet Lagrange-Funktion und 
% Lagrange-Gleichungen des Floating Arm Trebuchets 
% 
% Benutzt die symbolische Euler-Lagrange-Berechnung nach 
% Morten Veng (2021). Euler-Lagrange Solver 
% (https://www.mathworks.com/matlabcentral/fileexchange/93275-euler-lagrange-solver), 
% MATLAB Central File Exchange. Retrieved December 27, 2021.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Definiert die verallgemeinerten Koordinaten theta1, theta2 und Parameter
syms MG mw mT g delta gamma dgamma ddelta ls b B JB  'real'


%% Symbolische Lösung Lagrange-Funktion

% Verallgemeinerte Koordinaten und Ableitungen
q  = [gamma; delta];
dq = [dgamma; ddelta];

% kinetische Energien
xw = -(b+B)*cos(gamma) + ls*cos(delta);
zw = - b*sin(gamma)    + ls*sin(delta);
xG = 0;
zG = B*sin(gamma);

dxw = diff(xw,gamma)*dgamma + diff(xw,delta)*ddelta;
dzw = diff(zw,gamma)*dgamma + diff(zw,delta)*ddelta;
dxG = diff(xG,gamma)*dgamma + diff(xG,delta)*ddelta;
dzG = diff(zG,gamma)*dgamma + diff(zG,delta)*ddelta;

T_trans = 1/2* MG *(dxG.*dxG + dzG.*dzG) + 1/2 * JB *dgamma^2 + ...
          1/2* mw *(dxw.*dxw + dzw.*dzw);
T_rot = 0;
T_kin = T_trans + T_rot;

% potentielle Energien
U_pot =  g* MG * zG + mw * g * zw - 1/2 * mT * g * (b-B)*sin(gamma);

% Lagrange-Funktion
L = T_kin - U_pot;

%% Ausgabe der Ergebnisse
fprintf("------------------------------------\nLagrange-Funktion L:");
fprintf("\n------------------------------------\n");
fprintf(1, 'T = %s \n', T_kin);
fprintf(1, 'U = %s \n', U_pot);
fprintf(1, 'L = T - U = %s \n', L);
fprintf("\n");

%Externe Kräfte 
Q = [0; 0];

% Berechnet Euler-Lagrange-Gleichungen:
EQ = EulerLagrange(q,dq,L,Q,true);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
