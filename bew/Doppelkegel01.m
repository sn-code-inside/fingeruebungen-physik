% -------------------------------------------------------------------------
% Doppelkegel01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Aufwärtsrollender Doppelkegel
% 
% Programm berechnet Lagrange-Funktion und 
% Lagrange-Gleichungen des aufwärtsrollenden Doppelkegels  
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

syms M Js mw mT theta phi psi R H q dq alpha C1 C2 C3 g 'real'


%% Symbolische Lösung Lagrange-Funktion

% Verallgemeinerte Koordinaten und Ableitungen
Q  = [q];
dQ = [dq];

% kinetische Energie
% T_kin = 1/2* M *dq.*dq .* (1 + (3*R^2/10/tan(alpha)/tan(alpha))*...
%           ( 1/( H*tan(psi)/tan(alpha)-q ).^2));
T_kin = 1/2* dq.*dq .* (1 + C1./((C3-q).^2));
% potentielle Energien
% U_pot =  - M* g* q * sin(alpha-theta);
U_pot =  C2*q;
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
QF = [0];

% Berechnet Euler-Lagrange-Gleichungen:
EQ = EulerLagrange(Q,dQ,L,QF,true);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
