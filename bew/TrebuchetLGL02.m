% -------------------------------------------------------------------------
% TrebuchetLGL02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Trebuchet
%
% Programm berechnet Lagrange-Funktion und Lagrange-Gleichungen des
% Trebuchets 
% Symbolische Berechnung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Definiert die verallgemeinerten Koordinaten und Parameter
syms gamma dgamma beta dbeta delta ddelta h0 b B LG lS g ...
     JG JW JBG JbW MG mW J0 U0 JB lambda 'real'


%% Symbolische Lösung Lagrange-Funktion

% Verallgemeinerte Koordinaten und Ableitungen
q  = [gamma; beta; delta];
dq = [dgamma; dbeta; ddelta];

% kinetische Energien
xGB  = +B*cos(gamma);         % Aufhängung Gegengewicht
zGB  = +B*sin(gamma) ;
xG =  xGB + LG*cos(beta);     % Gegengewicht
zG =  zGB + LG*sin(beta);
xwB  = -b*cos(gamma);         % Aufhängung Wurfgewicht
zwB  = -b*sin(gamma);
xW   = xwB + lS*cos(delta);   % Wurfgewicht
zW   = zwB + lS*sin(delta);

% Differentiale der verallgemeinerten Koordinaten
dxW = diff(xW,gamma)*dgamma + diff(xW,delta)*ddelta;
dzW = diff(zW,gamma)*dgamma + diff(zW,delta)*ddelta;
dxG = diff(xG,gamma)*dgamma + diff(xG,beta)*dbeta;
dzG = diff(zG,gamma)*dgamma + diff(zG,beta)*dbeta;

% Kinetische Energien
T_kin = J0*dgamma*dgamma/2 +JG^2*dbeta^2/2 + JW^2*ddelta^2/2 +...
        JBG*dgamma*dbeta*cos(gamma-beta) - ...
        JbW*dgamma*ddelta*cos(gamma-delta);

% potentielle Energien
U_pot = U0*sin(gamma)+g*MG*LG*sin(beta)+g*mW*lS*sin(delta);

% ZB
ZB = h0-b*sin(gamma)+lS*sin(delta);

% Lagrange-Funktion
L = T_kin - U_pot +lambda*ZB;

fprintf("------------------------------------\nLagrange-Funktion L:");
fprintf("\n------------------------------------\n");
fprintf(1, 'T = %s \n', T_kin);
fprintf(1, 'U = %s \n', U_pot);
fprintf(1, 'L = T - U = %s \n', L);
fprintf("\n");

% Externe Kräfte 
Q = [0,0,0];

% Berechnet Euler-Lagrange-Gleichungen:
EQ = EulerLagrange(q,dq,L,Q,true);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


