% -------------------------------------------------------------------------
% SchaukelLG.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Lagrange Gleichungen Symbolische Berechnung
%
% Programm berechnet Lagrange-Funktion und Lagrange-Gleichungen der
% Schaukel 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Definiert die verallgemeinerten Koordinaten und Parameter

syms theta dtheta psi dpsi ...
     J1 J2 L MD mg g  'real'


%% Symbolische Lösung Lagrange-Funktion

% Verallgemeinerte Koordinaten und Ableitungen
q  = [theta; psi];
dq = [dtheta; dpsi];


% Kinetische Energien
T_kin = J1*dtheta^2/2+J2*(dtheta+dpsi)^2/2 - ...
        L*MD*dtheta*(dtheta+dpsi)*cos(psi);

% potentielle Energien
U_pot = -mg*L*g*cos(theta) + MD*g*cos(theta+psi);

% Lagrange-Funktion
L = T_kin - U_pot;

fprintf("------------------------------------\nLagrange-Funktion L:");
fprintf("\n------------------------------------\n");
fprintf(1, 'T = %s \n', T_kin);
fprintf(1, 'U = %s \n', U_pot);
fprintf(1, 'L = T - U = %s \n', L);
fprintf("\n");

% Externe Kräfte 
Q = [0,0];

% Berechnet Euler-Lagrange-Gleichungen:
EQ = EulerLagrange(q,dq,L,Q,true);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
