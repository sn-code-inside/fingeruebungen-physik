% -------------------------------------------------------------------------
% Doppelpendel01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Lagrange-Gleichung des Doppelpendels 
% 
% Benutzt die symbolische Euler/Lagrange Berechnung nach 
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

syms m1 m2 g L1 L2  th1 th2 dth1 dth2  'real'


%% Lagrange Funktion

% Verallgemeinerte Koordinaten und Ableitungen
q  = [th1; th2];
dq = [dth1; dth2];

% kinetische Energien
E_trans = 1/2*(m1+m2)*(L1*L1*dth1.*dth1) + ...
          1/2*m2*(L2*L2*dth2.*dth2)  + ...
          2*m1*L1*L2*dth1.*dth2.*cos(th1-th2);
E_rot = 0;

% potentielle Energien
E_pot = - (m1+m2)*g*L1*cos(th1) - m2*g*L2*cos(th2);

% Lagrange Funktion
L = E_trans + E_rot - E_pot;

% Externe Kräfte 
Q = [0; 0];

% Berechnet Euler-Lagrange-Gleichungen:
EQ = EulerLagrange(q,dq,L,Q,true);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


