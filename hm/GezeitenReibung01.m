% -------------------------------------------------------------------------
% GezeitenReibung01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Parameter der Gezeitenreibung auf der Erde
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
% alles in Einheiten von m und s
PlanetPara = ImportPlanetParameter('PlanetenParameter.csv')

oneday =  86400;
oneyear=  oneday*365.25;
mondsyn=  27.32*oneday;
TPE    =  oneday;
omegE  =  2*pi/TPE;
omegM  =  2*pi/mondsyn;
G      =  6.671e-11;
Masse  =  PlanetPara.Masse;

MMond  =  PlanetPara.Masse(12);
MErde  =  PlanetPara.Masse(3);
aEM    =  PlanetPara.aP(12);
RE     =  PlanetPara.RP(3);

k      =  0.2 ;         %Love-Zahl
delta  =  deg2rad(3);   %Vorlauf Gezeitenellipsoid

%%
% Berechnung
outstr(10,100)=' ';

% Slow Down Rate Erde Omega in 1/s
SlowDownE= -(45*omegM^2*MMond^2*RE^3*k*delta)/(4*omegE*MErde^2*aEM^3);
formatSpec = '\nSlow Down Rate Erde (Omega_E)     : %+5.2e%s';
outstr = sprintf(formatSpec, SlowDownE ,' 1/s ');
disp(outstr);
% Verlängerung Rotationsdauer Erde pro Jahrhundert in s
Delta_TE = -SlowDownE*TPE*oneyear*100;
formatSpec = 'Verlängerung Rotationsdauer Erde  : %+5.2e%s';
outstr = sprintf(formatSpec, Delta_TE ,' s pro Jahrhundert');
disp(outstr);
% Zeit bis zur synchronen Rotation in Jahren
tEnd     = abs(1/SlowDownE)/oneyear; 
formatSpec = 'Zeit bis zur synchronen Rotation  : %+5.2e%s\n';
outstr = sprintf(formatSpec, tEnd ,' Jahre');
disp(outstr);
% Bahnvergrößerungsrate Mond in s
Acc_aM   = -(4*MErde*RE^2*SlowDownE*omegE)/(5*MMond*aEM^2*omegM);
formatSpec = 'Bahnvergrößerungsrate Mond        : %+5.2e%s';
outstr = sprintf(formatSpec, Acc_aM ,' 1/s');
disp(outstr);
% Bahnvergrößerung Mond pro Jahrhundert in m
Delta_aM = Acc_aM*aEM*oneyear;
formatSpec = 'Bahnvergrößerungs Mond            : %+5.2e%s';
outstr = sprintf(formatSpec, Delta_aM ,' m pro Jahrhundert');
disp(outstr);
% Slow Down Rate Omega Mond in 1/s
SlowDownM= -3*Acc_aM/2;
formatSpec = 'Slow Down Rate (Omega M)          : %+5.2e%s';
outstr = sprintf(formatSpec, SlowDownM ,' 1/s');
disp(outstr);
% Verlängerung Synodischer Monat Mond pro Jahrhundert in s
Delta_TM = -SlowDownM*mondsyn*oneyear*100;
formatSpec = 'Verlängerung Synod. Monat Mond    : %+5.2e%s \n';
outstr = sprintf(formatSpec, Delta_TM ,' 1/s');
disp(outstr);

% Energiedissipation in W
AngMom   = (9*G*MMond^2*RE^6*delta*k)/(2*RE*aEM^6);
P        = AngMom*(omegM-omegE);
formatSpec = 'Energiedissipation                : %+5.2e%s \n';
outstr = sprintf(formatSpec, P ,' W');
disp(outstr);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------