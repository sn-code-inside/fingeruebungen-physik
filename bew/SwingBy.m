% -------------------------------------------------------------------------
% SwingBy.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Swing By Manöver am Jupiter
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%%
% Parameter

ME = 5.974e24;                      % Masse Erde in kg
MS = 1.988e30;                      % Masse Sonne in kg
MJ = 1.899e27;                      % Masse Jupiter in kg
G  = 6.67384e-11;                   % Gravitationskonstante m^3/kg/s^2
aE = 149*10^9;                      % Astronom. Einheit in m
aJ = 5.2.*aE;                       % Jupiter Bahnradius


vFL = sqrt(2*G*MS/aE);              % Fluchtgeschwindigkeit in m/s
VJ  = 2*pi*aJ/(11.9*365*86400)      % Bahngeschwindigkeit Jupiter in m/s
vR  = vFL *sqrt(1/5.2)              % Raumschiffgeschwindigkeit am Jupiter
vRJ = vFL/5.2                       % Projektin von vR auf Jupiterbahn
beta= acosd(vRJ/vR)                 % Winkel zwischen vR und VJ

uRs = sqrt(VJ^2+vR^2-2*VJ*vRJ)      % asymptotische Geschwindigkeit in Sigma'

uR  = uRs + VJ                      % asymptotische Geschwindigkeit in Sigma

