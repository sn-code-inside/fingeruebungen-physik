% -------------------------------------------------------------------------
% Traegheitsmomente.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Traegheitsmomente
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
% Variablen: Traegheitsmomente

%%Erde als Kugel und homogener Rotationsellipsoid

% Werte aus der Literatur
rhoE = 5.51*10^3;
JE   = 9.69e37;
ME   = 5.97e24;
RP   = 6356752;
RA   = 6378137;
R    = 6371000;

JK   = 2*4*pi*rhoE*R^5/15
JK   = 2*ME*R*R/5

JZ   = 2/5*ME*RA*RA
JX   = 1/5*ME*(RA*RA+RP*RP)

DelJ = JZ-JX
f1   = DelJ/JZ 
fE   = (RA-RP)/RA

% Zweischichtenmodell Erde (Kugel)

rhoM = 5.20*10^3;
rhoK = 11.0*10^3;
rK    = 2400000;  %innerer Radius, Kernradius


J    = 8*pi*(rhoK*rK^5 + rhoM*(R^5-rK^5))/15/JE;
M    = 4*pi/3*(rhoK*rK^3 + rhoM*(R^3-rK^3))/ME;

% ME2  = 4*pi*rhoE*R^3/3/ME


% Pendeluhr

mP = 0.016;
LP = 0.350;

JP = (1/3)*mP*LP^2*10^7

mKS= 0.156;
RKS= 0.055;
lKS= 0.257;

JKS= (1/2)*mKS*RKS^2+mKS*lKS^2*10^7

mR = 0.048;
LR = 0.140;
lR = 0.070;

JR = (1/12)*mR*LR^2 + mR*(lR+LR/2)^2*10^7


Jges = JR+JKS+JP

r=JKS/Jges
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------