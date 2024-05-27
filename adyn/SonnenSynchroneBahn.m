%%-------------------------------------------------------------------------
% SonnenSynchroneBahn.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Berechnung der Bahn und der Bodenspur eines sonnensynchronen Satelliten.
%
% -------------------------------------------------------------------------

%% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
GME = 398600.4415;  % in km^3/s^2
RE  = 6371;         % in km


%% VorwÃ¤rtsrechnung Orbit aus Bahnparametern

% Berechnung mit SSO Parametern
eP     = 0.00;
OmegaP = deg2rad(100.1213);
omegaP = deg2rad(0);
MP     = deg2rad(wrapTo360(0));  
K      =  3;  %Zyklus in Tage
N      = 43;  %Zyklus in Umrundungen der Erde
D_La1  = -360*K/N;  %geforderte Verschiebung aufsteigender Knoten
%Berechnet aus Bedingung N*(OmegaDot - ThetaDot)+TP = -K*360°
aP     = 7153.123;   %in km
J2     = 1.083e-3;  
iP     = deg2rad(98.498);
NPoints= 201;

%rate change of Omega per second 
OmegaDot  = -1.5*sqrt(GME/aP^3)*J2*cos(iP)*RE^2/aP^2;  % per second
OmegaDotd = OmegaDot*24*60*60;   %rate change of Omega per day 

%rate change of omega per day 
omegaDot  = -0.75*sqrt(GME/aP^3)*(1-5*(cos(iP))^2)*J2*RE^2/aP^2; % per second 
omegaDotd = omegaDot*24*60*60;   %rate change of Omega per day 

%rate change of mean anomaly
nDot  = -0.75*sqrt(GME/aP^3)*(1-3*(cos(iP))^2)*J2*RE^2/aP^2; % per second
nDotd = nDot*24*60*60;    %rate change of Omega per day
n0d   = sqrt(GME/aP^3)*86400;

%rate change of Theta
ThetaDotd = 360.9856473*2*pi/360; %rad per day
ThetaDot  = ThetaDotd/24/60/60;   %rad per second
MP        = wrapTo360(0);  % mean anomaly 

%Umlaufzeiten
TPn    = 2*pi*sqrt(aP^3/GME)/86400;  %Umlaufzeit Kepler-Bahn in d
TPkor  = 2*pi/(n0d+nDotd+omegaDotd); %geforderte Umlaufzeit für SSO

%Check for closure after N=3 and K=43
% x1 = N*(OmegaDotd-ThetaDotd)*2*pi/(n0d+nDotd+omegaDotd);
% x2 = -K*2*pi;
% Parameterübergabe

TSpan = 6/24; %Zeitspanne in h

% Berechnung Bahndaten zeitpunkt T1
dt1 = datetime('2020-01-01 00:00:00');
T1  = juliandate(dt1) - TSpan/24; % Julianisches Datum  ET
T1E = T1 + TSpan/24 + TPkor;
T_vector1 = linspace(T1,T1E,NPoints);  %in julianischen Tagen udn Bruchteilen
SatSSO.BaPa =[aP , eP, iP , MP, OmegaP, omegaP];
SatSSO.Name = 'SSO';
SatData1  = SatPQR_perturbed_orbit(T_vector1, SatSSO, GME, OmegaDotd, omegaDotd, nDotd);

% Berechnung Bahndaten zeitpunkt T2
T2  = T1 + 42*TPkor  - TSpan/24;  % Julianisches Datum  ET
T2E = T1 + 43*TPkor  + TSpan/24;
MP = MP + (n0d + nDotd + omegaDotd)*(T2-T1) ;
OmegaP = OmegaP+OmegaDotd*(T2-T1);
T_vector2 = linspace(T2,T2E,NPoints);  %in julianischen Tagen udn Bruchteilen
SatSSO.BaPa =[aP , eP, iP , MP, OmegaP, omegaP];
SatData2  = SatPQR_perturbed_orbit(T_vector2, SatSSO, GME, OmegaDotd, omegaDotd, nDotd);

% Umrechnung in geographische Breite und Länge
lat1(:) = rad2deg(SatData1.el);
theta1  = GMSTsat(T_vector1);
lon1(:) = wrapTo180(rad2deg(SatData1.az)-theta1);
lat2(:) = rad2deg(SatData2.el);
theta2  = GMSTsat(T_vector2);
lon2(:) = wrapTo180(rad2deg(SatData2.az)-theta2);


%% Graphik Ground Track

titlestr = strcat(SatSSO.Name,' Ground Track');
figure('name',titlestr);
gx = geoaxes;
hp(1)=geoplot(gx, lat1(:),lon1(:),'d',...
        'MarkerSize',2,'Color', Colors(2,:),'LineWidth', 3);
hold on
hp(2)=geoplot(gx, lat2(:),lon2(:),'+',...
        'MarkerSize',2,'Color', Colors(4,:),'LineWidth', 2);
geobasemap(gx,'bluegreen')
geolimits('manual') 
geolimits([-90 90],[-180 +180])
hp1 = title(titlestr,'FontSize',12);
legend(hp,'1. Orbit', '42. und 43. Orbit', 'location', 'northeast')
% legend box off
set(hp1,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',14);

%% PrintOut


fprintf('\n');
fprintf('Bahnparameter fÃ¼r Bahnberechnung nach Kepler\n');
fprintf('\n Zeit  = %s', dt1);
fprintf('\n');
fprintf('\n a     = %8.2f km  (GroÃŸe Halbachse)', aP);
fprintf('\n h     = %8.2f km  (mittlere Höhe)', aP-RE);
fprintf('\n e     = %8.5f     (ExzentrizitÃ¤t)', eP);
fprintf('\n i     = %8.2f     (Inklination)', iP);
fprintf('\n Omega = %8.2f     (Rektasz. aufst. Knoten)', rad2deg(OmegaP));
fprintf('\n omega = %8.2f     (Argument PerigÃ¤um)', rad2deg(omegaP));
fprintf('\n M     = %8.2f     (Mittlere Anomalie)', rad2deg(MP));
fprintf('\n TPn   = %8.2f min (Umlaufzeit Kepler-Bahn)', TPn*24*60);
fprintf('\n TPkor = %8.2f min (Umlaufzeit Sonnensynchrone Bahn)', TPkor*24*60);
fprintf('\n');
fprintf('\n D_La1 =%+8.2f min (Verschiebung aufsteigender Knoten gefordert)', D_La1);
fprintf('\n');
fprintf('\n');

%% Ende Programm


