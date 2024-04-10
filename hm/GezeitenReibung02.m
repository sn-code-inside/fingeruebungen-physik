% -------------------------------------------------------------------------
% GezeitenReibung02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Parameter der Gezeitenreibung auf der Erde
% nach dem Hantelmodell
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
rhoE   =  PlanetPara.rhoP(3);
k      =  0.20;         %Love-Zahl
delta  =  deg2rad(3);   %Vorlauf Gezeitenellipsoid

%%
% Berechnung
outstr(10,100)=' ';

% 1. Trägehitsmoment Erde
rhoK     =  13000;
rhoS     =   4000;
RK       =  3.5e06;
JE       =  8*pi*(rhoS*(RE^5-RK^5)+rhoK*RK^5)/15; %2 Schichten Modell
JEnaeh   =  0.2*MErde*RE^2;  %Näherung Kugel
ratio2   = JEnaeh/JE;

% 2. Drehimpulse  heute
JM       = MMond*aEM^2;
ratio1   = JEnaeh/JM;
LM0      = JM*omegM;
Lges     = LM0+ JE*omegE;

formatSpec = 'Gesamtdrehimpuls Erde+Mond-System : %+5.2e%s';
outstr = sprintf(formatSpec, Lges ,' Nm = const');
disp(outstr);
formatSpec = 'Trägheitsmoment Erde              : %+5.2e%s';
outstr = sprintf(formatSpec, JE ,' kg m^2 ');
disp(outstr);
formatSpec = 'Trägheitsmoment Mond              : %+5.2e%s';
outstr = sprintf(formatSpec, JM ,' kg m^2 ');
disp(outstr);

% 3. Berechnung der finalen Werte
aEMf     = Lges^2/(G*MErde*MMond^2);
omegMf   = Lges/MMond/aEMf^2;
ratio3   = aEMf/aEM;
TMf      = 2*pi/omegMf/86400;

formatSpec = 'Finale Bahnabstand Erde-Mond      : %+5.2e%s';
outstr = sprintf(formatSpec, aEMf ,' m ');
formatSpec = 'Finale Rotationsdauer Erde = Mond : %+5.2f%s';
outstr = sprintf(formatSpec, TMf ,' d ');
disp(outstr);

LMf      = MMond*aEMf^2*omegMf;
LEf      = JE*omegMf;
ratio4   = LEf/LMf;

% 4. Drehmoment
m        = MErde*1e-08; %Hantelmasse
delta    = deg2rad(3);  %Verdrehung des Gezeitenellipsoids
FA       = G*m*MMond/(aEM^2+RE^2-2*aEM*RE*cos(delta));
FB       = G*m*MMond/(aEM^2+RE^2+2*aEM*RE*cos(delta));
MDges    = RE*delta*(FB-FA); %Drehmoment
formatSpec = 'Gesamtdrehmoment auf Hantel       : %+5.2e%s';
outstr = sprintf(formatSpec, MDges ,' Nm ');
disp(outstr);

d_omegE  = MDges/JE;
ratio5   = d_omegE/omegE;
d_aEM    = -2*JE*d_omegE/(MMond*aEM*omegM);
ratio6   = -d_aEM/aEM;

% 5. Bahnvergrößerungsrate Mond in s
DeltaaEM = d_aEM*100*365*86400;
formatSpec = 'Bahnänderungssrate(d_aEM/d_t)/aEM : %+5.2e%s';
outstr = sprintf(formatSpec, d_aEM/aEM ,' 1/s');
disp(outstr);
% Bahnvergrößerung Mond pro Jahrhundert in m
formatSpec = 'Vergrößerung Mondbahn             : %+5.2f%s';
outstr = sprintf(formatSpec, DeltaaEM ,' m pro Jahrhundert');
disp(outstr);
% 6. Verlängerung Rotationsdauer Erde
DeltaTPE = -(d_omegE*2*pi*100*365*86400)/omegE^2;
formatSpec = 'Verlängerung Rotationsdauer Erde  : %+5.2e%s';
outstr = sprintf(formatSpec, DeltaTPE ,' s pro Jahrhundert');
disp(outstr);

% 6. Energiebetrachtungen
Eges     = JE*omegE^2/2*JM*omegM^2/2-G*MMond*MErde/aEM;
Pverl    = JE*omegE*d_omegE-1.5*JM*G*MErde*d_aEM/aEM^4 ...
           + G*MErde*MMond*d_aEM/aEM^2;
formatSpec = '\nGesamtenergie System              : %+5.2e%s';
outstr = sprintf(formatSpec, Eges ,' Ws');
disp(outstr);
formatSpec = 'Verlustleistung                   : %+5.2e%s';
outstr = sprintf(formatSpec, Pverl ,' W');
disp(outstr);

% 7. Energiedissipation in Wassersäulenmodell
hw       = 0.5;  %Höhe Wassersäule
rhow     = 1000; %Dichte Wassersäule
g        = 9.81;     
MW       = (4/3)*pi*rhow*RE^3*3*hw/RE; % Masse      
EWperDay = -2*MW*g*hw; %Arbeit durch Gezeiten
PW       = 0.1*EWperDay/86400;  %Leistung
EWperJ   = PW*365*86400/3600;  %in kWh
ratio7   = PW/Pverl;
formatSpec = 'Verlustleistung (Wassersäulem.)   : %+5.2e%s';
outstr = sprintf(formatSpec, PW ,' W');
disp(outstr);
formatSpec = 'Verlustleistung (Wassersäulem.)   : %+8.3f%s';
outstr = sprintf(formatSpec, PW*1e-09 ,' GW');
disp(outstr);
formatSpec = 'Energieverlust p.a.(Wassersäulem.): %+5.2e%s';
outstr = sprintf(formatSpec, EWperJ ,' kWh');
disp(outstr);
formatSpec = 'Energieverlust p.a.(Wassersäulem.): %+5.2e%s';
outstr = sprintf(formatSpec, EWperJ*1e-09 ,' TWh');
disp(outstr);

% 8. Erwärmung Ozeane pro Jahr  
cW       = 4.2e3;  %spezifische Wärmekapazität
DeltaTW  = -PW/MW/cW;       
DeltaTWpJ= -EWperDay*365/MW/cW;
formatSpec = 'Erwärmung Ozeane pro Jahr         : %+5.2f%s';
outstr = sprintf(formatSpec, DeltaTWpJ ,' K');
disp(outstr);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------